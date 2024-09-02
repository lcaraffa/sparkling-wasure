# Convertcloud

Simple pointcloud format converter. Supported conversions:


|format|read|write|ascii|binary|compressed|
|------|----|------|----|------|----------|
| .pcd |![green]|![green]|![green]|![red]|![red]|
| .ply |![green]|![green]|![green]|![red]|![gray]|
| .xyz |![green]|![green]|![green]|![gray]|![gray]|
| .zdf |![green]|![red]|![gray]|![gray]|![green]|
| .a3d |![green]|![green]|![green]|![gray]|![gray]|
| .stl |![green]|![red]|![green]|![red]|![gray]|

Use from shell: 
```sh
$ cvc original.format1 converted.format2 
```

Use from python:
```python
import convertcloud as cvc

conv = cvc.Converter()
conv.load_points("original.format1")
conv.convert("converted.format2")
```
[green]: https://via.placeholder.com/20/36b023/?text=+
[red]: https://via.placeholder.com/20/f03c15/?text=+
[gray]: https://via.placeholder.com/20/c4c4c4/?text=+
