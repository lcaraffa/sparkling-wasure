python -m pdaltools.replace_attribute_in_las \
    --input_file test/data/classified_laz/test_data_77050_627755_LA93_IGN69.laz \
    --output_file test/tmp/replaced_cmdline.laz \
    --replacement_map test/data/example_replacement_map.json \
    --record_format 8