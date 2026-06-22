import json
import os
import struct
import sys

GLB_MAGIC = b"glTF"
GLB_CHUNK_JSON = 0x4E4F534A
GLB_CHUNK_BIN = 0x004E4942


def rewrite_contents(node):
    """Walk the global tileset and turn each leaf into an *external tileset*
    reference: ``content.uri`` is rewritten from ``tiles/<hash>.obj`` to
    ``<hash>/tileset.json`` (the autonomous, ECEF-georeferenced tileset that
    py3dtilers produced for that tile). Also strips two non-spec artefacts
    written by mesh23dtile.py: the custom ``raw_bbox`` field and empty
    ``children`` arrays (3D Tiles 1.0 forbids ``children: []``).
    """
    if isinstance(node, dict):
        # Drop the non-standard helper field and any empty children array.
        node.pop("raw_bbox", None)
        if isinstance(node.get("children"), list) and len(node["children"]) == 0:
            del node["children"]
        # Redirect the leaf URI to the external tileset.
        content = node.get("content")
        if isinstance(content, dict) and "uri" in content:
            bname = os.path.basename(content["uri"])
            fname, _ = os.path.splitext(bname)
            content["uri"] = fname + "/tileset.json"
        for value in node.values():
            rewrite_contents(value)
    elif isinstance(node, list):
        for item in node:
            rewrite_contents(item)


def fix_b3dm_alpha(path):
    """Force ``baseColorFactor[3] = 1.0`` on OPAQUE materials whose alpha is
    written as 0 by py3dtilers v9 (otherwise the mesh is invisible in iTowns
    even though ``alphaMode = "OPAQUE"``).

    The function looks long because it edits a binary file embedded inside
    another binary file:

        b3dm  ::=  b3dm header (28 B)
                || feature table JSON || feature table BIN
                || batch table JSON   || batch table BIN
                || GLB
        GLB   ::=  GLB header (12 B)
                || JSON chunk header (8 B) || JSON chunk data
                || BIN chunk header  (8 B) || BIN chunk data

    Changing one field inside the glTF JSON shifts all the following bytes, 
    so we have to rebuild the GLB and then the b3dm.
    """
    # --- 1. Read the file and sanity-check the b3dm magic.
    raw = open(path, "rb").read()
    if len(raw) < 28 or raw[:4] != b"b3dm":
        return

    # --- 2. Parse the 28-byte b3dm header to find each section length.
    (
        _magic,
        _version,
        _byte_length,
        ft_json_len,
        ft_bin_len,
        bt_json_len,
        bt_bin_len,
    ) = struct.unpack_from("<4sIIIIII", raw, 0)

    # --- 3. Split the file into "everything before the GLB" (kept verbatim)
    # and the GLB itself (the only part we will touch).
    pre_glb_size = 28 + ft_json_len + ft_bin_len + bt_json_len + bt_bin_len
    pre_glb = raw[:pre_glb_size]
    glb = raw[pre_glb_size:]
    if glb[:4] != GLB_MAGIC:
        return

    # --- 4. Parse the GLB header (12 B) and decode the JSON chunk.
    # Strip both legal paddings (NULs *and* spaces) before json.loads.
    glb_version, _glb_length = struct.unpack_from("<II", glb, 4)
    json_chunk_len, json_chunk_type = struct.unpack_from("<II", glb, 12)
    if json_chunk_type != GLB_CHUNK_JSON:
        return
    json_bytes = glb[20 : 20 + json_chunk_len]
    gltf = json.loads(json_bytes.decode("utf-8").rstrip("\x00 "))

    # --- 5. Extract the BIN chunk. We do NOT modify it, but we need to
    # carry it over to the rebuilt GLB and may have to extend it later for
    # 8-byte alignment.
    bin_off = 20 + json_chunk_len
    bin_chunk_len, bin_chunk_type = struct.unpack_from("<II", glb, bin_off)
    if bin_chunk_type != GLB_CHUNK_BIN:
        return
    bin_data = glb[bin_off + 8 : bin_off + 8 + bin_chunk_len]

    # --- 6. The actual semantic patch: any OPAQUE material whose
    # baseColorFactor has alpha=0 gets its alpha set to 1. If nothing
    # matched, exit early and leave the file untouched.
    edited = False
    for mat in gltf.get("materials", []):
        if mat.get("alphaMode", "OPAQUE") != "OPAQUE":
            continue
        pbr = mat.get("pbrMetallicRoughness")
        if not pbr:
            continue
        bcf = pbr.get("baseColorFactor")
        if isinstance(bcf, list) and len(bcf) == 4 and bcf[3] == 0:
            bcf[3] = 1.0
            edited = True

    if not edited:
        return

    # --- 7. Re-encode the JSON chunk.
    new_json = json.dumps(gltf, separators=(",", ":")).encode("utf-8")
    new_json += b" " * ((4 - (len(new_json) % 4)) % 4)
    new_bin = bin_data

    # --- 8. Rebuild the GLB from scratch with the new JSON chunk and the
    # untouched BIN chunk. The GLB total length is 12 (header) + 8 (JSON
    # chunk header) + JSON + 8 (BIN chunk header) + BIN.
    new_glb_length = 12 + 8 + len(new_json) + 8 + len(new_bin)
    new_glb = bytearray()
    new_glb += GLB_MAGIC
    new_glb += struct.pack("<II", glb_version, new_glb_length)
    new_glb += struct.pack("<II", len(new_json), GLB_CHUNK_JSON)
    new_glb += new_json
    new_glb += struct.pack("<II", len(new_bin), GLB_CHUNK_BIN)
    new_glb += new_bin

    # --- 9. Enforce the 3D Tiles rule that the whole b3dm must be 8-byte aligned.
    new_total = pre_glb_size + len(new_glb)
    while new_total % 8 != 0:
        extra = 8 - (new_total % 8)
        new_bin = new_bin + b"\x00" * extra
        gltf["buffers"][0]["byteLength"] = len(new_bin)
        new_json = json.dumps(gltf, separators=(",", ":")).encode("utf-8")
        new_json += b" " * ((4 - (len(new_json) % 4)) % 4)
        new_glb_length = 12 + 8 + len(new_json) + 8 + len(new_bin)
        new_glb = bytearray()
        new_glb += GLB_MAGIC
        new_glb += struct.pack("<II", glb_version, new_glb_length)
        new_glb += struct.pack("<II", len(new_json), GLB_CHUNK_JSON)
        new_glb += new_json
        new_glb += struct.pack("<II", len(new_bin), GLB_CHUNK_BIN)
        new_glb += new_bin
        new_total = pre_glb_size + len(new_glb)

    # --- 10. Patch the b3dm header's byteLength (offset 8, 4 bytes LE) so
    # it matches the new total size, then write the file back.
    new_header = bytearray(pre_glb)
    struct.pack_into("<I", new_header, 8, new_total)
    with open(path, "wb") as f:
        f.write(bytes(new_header) + bytes(new_glb))


def patch_all_b3dm(tile_path):
    """Applies fix_b3dm_alpha to all .b3dm files found under tile_path."""
    for dirpath, _dirnames, filenames in os.walk(tile_path):
        for name in filenames:
            if name.endswith(".b3dm"):
                fix_b3dm_alpha(os.path.join(dirpath, name))


def main(tile_path):
    in_path = os.path.join(tile_path, "tileset_tmp.json")
    out_path = os.path.join(tile_path, "tileset.json")
    with open(in_path, "r") as file:
        data = json.load(file)
    rewrite_contents(data)
    with open(out_path, "w") as file:
        json.dump(data, file, indent=4)
    patch_all_b3dm(tile_path)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python finalize.py <output_dir>")
        sys.exit(1)
    main(sys.argv[1])
