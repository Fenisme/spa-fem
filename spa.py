#!/usr/bin/env freecadcmd
import FreeCAD as App
from FreeCAD import Vector, Rotation
from Part import LineSegment
from Sketcher import Constraint
from PartDesign import *
from ObjectsFem import *
from Fem import FemMesh
import shutil

class Spa:
    left = 0
    right = 1
    up = 2
    down = 3
    horizontal = 0
    vertical = 1
    start = 1
    end = 2
    # n: number of chamber
    # w: width
    # h: height
    # a, dn: chamber length, length between chamber
    # t, tb, tc: thickness of chamber, bottom, air channel
    # lp, le: length of start, end
    def __init__(self, n: int, w: float, h: float, a: float, dn: float, t: float, tb: float, tc: float, lp: float, le: float, tu: float = 0):
        self.n = n
        self.w = w
        self.h = h
        self.a = a
        self.dn = dn
        self.t = t
        self.tb = tb
        self.tc = tc
        self.lp = lp
        self.le = le
        if tu == 0:
            self.tu = self.t
        else:
            self.tu = tu

    def GetDoc(self):
        if App.activeDocument() is None:
            App.newDocument()
        return App.activeDocument()

    def __lineGetPoint(self, sketch, index: int, direction: int):
        line = sketch.Geometry[index]
        if direction == Spa.up:
            if line.StartPoint.y > line.EndPoint.y:
                return self.start
            elif line.StartPoint.y < line.EndPoint.y:
                return self.end
        elif direction == Spa.down:
            if line.StartPoint.y > line.EndPoint.y:
                return self.end
            elif line.StartPoint.y < line.EndPoint.y:
                 return self.start
        elif direction == Spa.right:
            if line.StartPoint.x > line.EndPoint.x:
                return self.start
            elif line.StartPoint.x < line.EndPoint.x:
                return self.end
        elif direction == Spa.left:
            if line.StartPoint.x > line.EndPoint.x:
                return self.end
            elif line.StartPoint.x < line.EndPoint.x:
                return self.start
        print("Warning: direction is seem to be wrong")
        return self.start

    def __lineGetPointVector(self, sketch, index: int, direction: int):
        if self.__lineGetPoint(sketch, index, direction) == Spa.start:
            return sketch.Geometry[index].StartPoint
        else:
            return sketch.Geometry[index].EndPoint

    def __sketcherAddLine(self, sketch, length: float, src_index: int, src_point: int, direction: int):
        line_index = sketch.addGeometry(LineSegment(Vector(0, 0, 0), Vector(1, 1, 0)), False)
        sketch.addConstraint(Constraint("Distance", line_index, length))
        if direction in (Spa.up, Spa.down):
            sketch.addConstraint(Constraint("Vertical", line_index))
        elif direction in (Spa.left, Spa.right):
            sketch.addConstraint(Constraint("Horizontal", line_index))
        else:
            print("__sketcherAddLine receive an unknown direction:", direction)
            return

        sketch.movePoint(line_index, self.__lineGetPoint(sketch, line_index, direction), self.__lineGetPointVector(sketch, src_index, src_point))
        sketch.addConstraint(Constraint("Coincident", line_index, self.__lineGetPoint(sketch, line_index, direction), src_index, self.__lineGetPoint(sketch, src_index, src_point)))
        return line_index

    def __getFace(self, solid, point):
        index = 0
        for i in range(len(solid.Faces)):
            if solid.Faces[i].Surface.projectPoint(point, "LowerDistance") < solid.Faces[index].Surface.projectPoint(point, "LowerDistance"):
                index = i
            if solid.Faces[i].Surface.projectPoint(point, "LowerDistance") == 0:
                return i
        print("__getFace: Cannot find a face match the point exactly")
        return index

    def __getFaces(self, solid, point):
        lst = []
        for i in range(len(solid.Faces)):
             if solid.Faces[i].Surface.projectPoint(point, "LowerDistance") < 0.001:
                 lst.append(i)
        return lst

    def __getFaceName(self, solid, point):
        return "Face%d"%(self.__getFace(solid, point) + 1)

    def __addFaces(self, lst: list, face: list):
        for i in face:
            if i not in lst:
                lst.append(i)

    def makeOuterSketcher(self, sketch):
        #sketch = self.GetDoc().addObject("Sketcher::SketchObject", name)
        sketch.Placement = App.Placement(Vector(0, 0, 0), Rotation(0, 0, 0, 1))

        # 0: Bottom line
        index_last = sketch.addGeometry(LineSegment(Vector(0, 0, 0), Vector(1, 0, 0)), False)
        sketch.addConstraint(Constraint("Horizontal", 0))
        sketch.addConstraint(Constraint("Distance", 0, self.lp + self.le + self.n * self.a + (self.n - 1) * self.dn))
        sketch.addConstraint(Constraint("Coincident", 0, self.__lineGetPoint(sketch, 0, self.left), -1, 1))

        # Loop struct is |_|-, deal with lp first
        if self.lp == 0:
            index_last = self.__sketcherAddLine(sketch, self.h, index_last, self.left, self.down)
        else:
            index_last = self.__sketcherAddLine(sketch, self.tb + self.tc + self.t, index_last, self.left, self.down)
            index_last = self.__sketcherAddLine(sketch, self.lp, index_last, self.up, self.left)
            index_last = self.__sketcherAddLine(sketch, self.h - self.tb - self.tc - self.t, index_last, self.right, self.down)

        index_last = self.__sketcherAddLine(sketch, self.a, index_last, self.up, self.left)

        for i in range(1, self.n):
            index_last = self.__sketcherAddLine(sketch, self.h - self.tb - self.tc - self.t, index_last, Spa.right, Spa.up)
            index_last = self.__sketcherAddLine(sketch, self.dn, index_last, Spa.down, Spa.left)
            index_last = self.__sketcherAddLine(sketch, self.h - self.tb - self.tc - self.t, index_last, Spa.right, Spa.down)
            index_last = self.__sketcherAddLine(sketch, self.a, index_last, Spa.up, Spa.left)

        if self.le != 0:
            index_last = self.__sketcherAddLine(sketch, self.h - self.tb - self.tc - self.t, index_last, Spa.right, Spa.up)
            index_last = self.__sketcherAddLine(sketch, self.le, index_last, Spa.down, Spa.left)

        # last: connect two point
        index_last = sketch.addGeometry(LineSegment(Vector(0, 0, 0), Vector(0, 1, 0)), False)
        sketch.addConstraint(Constraint("Coincident", 0, self.__lineGetPoint(sketch, 0, self.right), index_last, self.start))
        sketch.addConstraint(Constraint("Coincident", index_last - 1, self.__lineGetPoint(sketch, index_last - 1, self.right), index_last, self.end))

        return sketch

    def __makeAirchannelSketcher(self, sketch):
        # 0: Bottom line
        index_last = sketch.addGeometry(LineSegment(Vector(-(self.tb + self.tc), (self.w - self.tc) / 2, 0), Vector(-self.tb, (self.w - self.tc) / 2, 0)), False)
        sketch.addConstraint(Constraint("Horizontal", index_last))
        sketch.addConstraint(Constraint("Distance", index_last, self.tc))
        sketch.addConstraint(Constraint("DistanceX", 0, self.__lineGetPoint(sketch, 0, self.right), -1, 1, self.tb))
        sketch.addConstraint(Constraint("DistanceY", -1 , 1, index_last, self.__lineGetPoint(sketch, 0, self.right), (self.w - self.tc) / 2))

        index_last = self.__sketcherAddLine(sketch, self.tc, index_last, Spa.left, Spa.down)
        index_last = self.__sketcherAddLine(sketch, self.tc, index_last, Spa.up, Spa.left)

        # last: connect two point
        index_last = sketch.addGeometry(LineSegment(Vector(0, 0, 0), Vector(0, 1, 0)), False)
        sketch.addConstraint(Constraint("Coincident", 0, self.__lineGetPoint(sketch, 0, self.right), index_last, self.start))
        sketch.addConstraint(Constraint("Coincident", index_last - 1, self.__lineGetPoint(sketch, index_last - 1, self.right), index_last, self.end))

    def __makeChamberSketcher(self, sketch):
        for i in range(self.n):
            # bottom
            index_last = sketch.addGeometry(LineSegment(Vector(-(self.lp + self.t + i * (self.a + self.dn)), self.t, 0), Vector(-(self.lp + self.a - self.t + i * (self.a + self.dn)), self.t, 0)), False)
            sketch.addConstraint(Constraint("Horizontal", index_last))
            sketch.addConstraint(Constraint("Distance", index_last, self.a - 2 * self.t))
            if i == 0:
                sketch.addConstraint(Constraint("DistanceX", index_last, self.__lineGetPoint(sketch, index_last, self.right), -1, 1, self.lp + self.t + i * (self.a + self.dn)))
                sketch.addConstraint(Constraint("DistanceY", -1, 1, index_last, self.__lineGetPoint(sketch, index_last, self.right), self.t))
            else:
                sketch.addConstraint(Constraint("DistanceX", index_last, self.__lineGetPoint(sketch, index_last, self.right), index_last - 4, self.__lineGetPoint(sketch, index_last - 4, self.right), self.a + self.dn))
                sketch.addConstraint(Constraint("DistanceY", index_last, self.__lineGetPoint(sketch, index_last, self.right), index_last - 4, self.__lineGetPoint(sketch, index_last - 4, self.right), 0))
            
            index_last = self.__sketcherAddLine(sketch, self.w - 2 * self.t, index_last, self.left, self.down)
            index_last = self.__sketcherAddLine(sketch, self.a - 2 * self.t, index_last, self.up, self.left)

            index_last = sketch.addGeometry(LineSegment(Vector(0, 0, 0), Vector(0, 1, 0)), False)
            sketch.addConstraint(Constraint("Coincident", index_last - 1, self.__lineGetPoint(sketch, index_last - 1, self.right), index_last, self.start))
            sketch.addConstraint(Constraint("Coincident", index_last - 3, self.__lineGetPoint(sketch, index_last - 3, self.right), index_last, self.end))

    def makeBody(self, name:str = 'Body_spa'):
        body = self.GetDoc().addObject("PartDesign::Body", name)
        
        # Outer shape
        sketch_outer = body.newObject("Sketcher::SketchObject", "Skecher_outer")
        sketch_outer.Support = (self.GetDoc().getObject('XY_Plane'), [''])
        sketch_outer.MapMode = "FlatFace"
        self.makeOuterSketcher(sketch_outer)
        pad = body.newObject("PartDesign::Pad", "Pad_spa")
        pad.Profile = sketch_outer
        pad.Length = self.w
        pad.ReferenceAxis = (sketch_outer, ['N_Axis'])
        self.GetDoc().recompute()

        # Air channel
        sketch_inner = body.newObject("Sketcher::SketchObject", "Sketch_air_channel")
        sketch_inner.Support = (pad, [self.__getFaceName(pad.Shape, Vector(0, (self.tb + self.tc + self.t) / 2, self.w / 2))])
        sketch_inner.MapMode = "FlatFace"
        self.__makeAirchannelSketcher(sketch_inner)
        pocket = body.newObject("PartDesign::Pocket", "Pocket_spa")
        pocket.Profile = sketch_inner
        pocket.Length = self.lp + self.n * self.a + (self.n - 1) * self.dn + self.le - self.t
        pocket.ReferenceAxis = (sketch_inner, ['N_Axis'])
        self.GetDoc().recompute()

        # Chambers
        sketch_chambers = body.newObject("Sketcher::SketchObject", "Sketch_chambers")
        sketch_chambers.Support = (pocket, [self.__getFaceName(pocket.Shape, Vector((self.lp + self.n * self.a + (self.n - 1) * self.dn + self.le - self.t) / 2, self.tb, self.w / 2))])
        sketch_chambers.MapMode = "FlatFace"
        self.__makeChamberSketcher(sketch_chambers)
        pocket = body.newObject("PartDesign::Pocket", "Pocket_chambers")
        pocket.Profile = sketch_chambers
        pocket.Length = self.h - self.tb - self.t - self.tu
        pocket.ReferenceAxis = (sketch_chambers, ['N_Axis'])
        pocket.Reversed = 1
        self.GetDoc().recompute()

        return body

    def makeMaterial(self):
        if not self.ccx_working_dir:
            print(".inp file not found")
            return
        shutil.copy(self.ccx_working_dir + "/Mesh.inp", self.ccx_working_dir + "/Mesh-old.inp")
        f_in = open(self.ccx_working_dir + "/Mesh-old.inp", 'r')
        if not f_in:
            print(".inp file cannot open")
            return
        f_out = open(self.ccx_working_dir + "/Mesh.inp", 'w')
        if not f_out:
            print(".inp file cannot write")
            return

        n_del = 0
        for i in f_in.readlines():
            if "*ELASTIC" in i:
                print("*HYPERELASTIC,YEOH", file = f_out)
                #print("573.82,-74.744,-11.321,0.01,0.1,0.5", file = f_out)
                print("110,20,0,0,0,0", file = f_out)
                n_del = 1
                continue
            elif "*STEP" in i:
                print("*STEP, NLGEOM", file = f_out)
                continue
            elif n_del > 0:
                n_del -= 1
                continue
            print(i, end = '', file = f_out)
        f_in.close()
        f_out.close()

    def makeSpaAnalysis(self, pressure:float = 0.5):
        spa = self.makeBody()
        # Analysis
        analysis = makeAnalysis(self.GetDoc(), "Analysis_spa")

        # Material
        material_obj = makeMaterialSolid(self.GetDoc(), "MechanicalMaterial")
        mat = material_obj.Material
        mat["Name"] = "Spa Silica"
        mat["YoungsModulus"] = "200000 MPa"
        mat["PoissonRatio"] = "0.30"
        mat["Density"] = "2200 kg/m^3"
        material_obj.Material = mat
        analysis.addObject(material_obj)

        # Constranint fixed
        con_fixed = makeConstraintFixed(self.GetDoc(), "ConstranintFixed")
        con_fixed.References = [(spa, self.__getFaceName(spa.Shape, Vector(0, self.tb / 2, self.w / 2)))]
        analysis.addObject(con_fixed)

        # Pressure constranint
        con_pressure = makeConstraintPressure(self.GetDoc(), "ConstranintPressure")
        refs_index = []
        # Air channel
        self.__addFaces(refs_index, self.__getFaces(spa.Shape, Vector(self.lp + self.n * self.a + (self.n - 1) * self.dn + self.le - self.t, self.tb + self.tc / 2, self.w / 2)))
        self.__addFaces(refs_index, self.__getFaces(spa.Shape, Vector((self.lp + self.t) / 2, self.tb, self.w / 2)))
        self.__addFaces(refs_index, self.__getFaces(spa.Shape, Vector((self.lp + self.t) / 2, self.tb + self.tc, self.w / 2)))
        self.__addFaces(refs_index, self.__getFaces(spa.Shape, Vector((self.lp + self.t) / 2, self.tb + self.tc / 2, self.w / 2 - self.tc / 2)))
        self.__addFaces(refs_index, self.__getFaces(spa.Shape, Vector((self.lp + self.t) / 2, self.tb + self.tc / 2, self.w / 2 + self.tc / 2)))
        
        # Chambers
        self.__addFaces(refs_index, self.__getFaces(spa.Shape, Vector(self.lp + self.t + self.a / 2, self.h - self.t, self.w / 2)))
        for i in range(self.n):
            self.__addFaces(refs_index, self.__getFaces(spa.Shape, Vector(self.lp + self.t + i * (self.a + self.dn), (self.h - self.t + self.tb + self.tc) / 2, self.w / 2)))
            self.__addFaces(refs_index, self.__getFaces(spa.Shape, Vector(self.lp + self.t + self.a - 2 * self.t + i * (self.a + self.dn), (self.h - self.t + self.tb + self.tc) / 2, self.w / 2)))
            self.__addFaces(refs_index, self.__getFaces(spa.Shape, Vector(self.lp + self.a / 2, (self.h - self.t + self.tb + self.tc) / 2, self.t)))
            self.__addFaces(refs_index, self.__getFaces(spa.Shape, Vector(self.lp + self.a / 2, (self.h - self.t + self.tb + self.tc) / 2, self.w - self.t)))
        refs = []
        for i in refs_index:
            refs.append((spa, "Face%d"%(i + 1)))
        con_pressure.References = refs
        con_pressure.Pressure = pressure
        analysis.addObject(con_pressure)

        # Mesh
        femmesh_obj = analysis.addObject(makeMeshGmsh(self.GetDoc(), "Mesh"))[0]
        femmesh_obj.Part = spa
        femmesh_obj.CharacteristicLengthMax = 10
        from femmesh import gmshtools
        gmsh_mesh = gmshtools.GmshTools(femmesh_obj, analysis)
        gmsh_mesh.create_mesh()
        self.gmsh_working_dir = gmsh_mesh.working_dir
        self.GetDoc().recompute()

        # Solver
        solver_obj = makeSolverCalculixCcxTools(self.GetDoc(), "CalculiXCcxTools")
        analysis.addObject(solver_obj)
        from femtools import ccxtools
        fea = ccxtools.FemToolsCcx()
        fea.update_objects()
        fea.setup_working_dir()
        self.ccx_working_dir = fea.working_dir
        fea.setup_ccx()
        message = fea.check_prerequisites()
        if not message:
            fea.purge_results()
            fea.write_inp_file()
            self.makeMaterial()
            fea.ccx_run()
            fea.load_results()
        else:
            print(message)

        return analysis

    def getAngle(self):
        try:
            results = self.GetDoc().getObject("CCX_Time_1_0_Results")
        except:
            results = self.GetDoc().getObject("CCX_Results")
        print("Min y:", results.Stats[2])
        return results.Stats[2]

    def clean(self):
        print("Cleaning", self.gmsh_working_dir, self.ccx_working_dir)
        if "tmp" in self.gmsh_working_dir:
            shutil.rmtree(self.gmsh_working_dir)
        if "tmp" in self.ccx_working_dir:
            shutil.rmtree(self.ccx_working_dir)
        App.closeDocument("Unnamed")

def makeData():
    la = []
    lb = []
    f = open("data.csv", 'r')
    for i in f.readlines():
        if ',' in i:
            a, b = eval(i.strip('\n'))
            la.append(a)
            lb.append(b)
    f.close()
    f = open("data.csv", 'w')
    for i in range(len(la)):
        print("{},{}".format(la[i],lb[i]), file = f)


    for i in range(200, 400):
        if i in la:
            continue
        spa = Spa(n = 6, w = i / 10, h = 3 + 20 + 1.5 + 4, a = 10, dn = 4, t = 1.5, tb = 10, tc = 2, lp = 0, le = 10)
        spa.makeSpaAnalysis()
        print("{},{}".format(i, spa.getAngle()), file = f)
        print("{},{}".format(i, spa.getAngle()))
        spa.clean()
        del(spa)
    f.close()

makeData()
#spa = Spa(n = 6, w = 30, h = 3 + 20 + 1.5 + 4, a = 10, dn = 4, t = 1.5, tb = 10, tc = 2, lp = 4, le = 10)
#spa.makeSpaAnalysis()
#spa.getAngle()
#p = spa.makeBody()

