# FEM for Soft Pneumatic Actuator
## Requirements

- FreeCAD
- Gmsh
- CalculiX-ccx

The script is only tested on Linux.
## Usage
If you want to run in terminal, use ```freecadcmd```:

```
freecadcmd spa.py
```

If you want to run in FreeCAD and get 3D modeling, run in Python console:

```
exec(open("spa.py", 'r').read())
```

You should edit spa.py manually to do 3D modeling or FEM analysis.

If you want to do 3D modeling, comment following line and create a Spa class, then call the ```make_body```

```
#makeData()
```

If you want to do FEM analysis, call ```makeSpaAnalysis```, or edit body of ```makeData``` to get a ```.csv``` result.

