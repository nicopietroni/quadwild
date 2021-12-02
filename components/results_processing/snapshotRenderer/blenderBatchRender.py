import bpy
import glob, os, sys

argv = sys.argv
argv = argv[argv.index("--") + 1:]

directory = os.path.abspath(argv[0])
patternglob = argv[1]
imageDir  = os.path.abspath(argv[2])

print(directory)
print(imageDir)

pattern = os.path.join(directory, patternglob)

print (pattern)

scene = bpy.context.scene
data  = bpy.data

for file in glob.glob(pattern):
    print(file)
    
    name = os.path.basename(file)
    
    success = bpy.ops.import_scene.obj(filepath=file, axis_forward='-Z', axis_up='Y')

    # the imported objects...
    obj_objects = bpy.context.selected_objects[:]

    if (len(obj_objects) == 0):
        continue

    obj = bpy.context.selected_objects[0]

    bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS', center='MEDIAN')
    obj.location = (0, 0, 0)
    obj.scale = obj.scale / max(obj.dimensions)

    for edge in obj.data.edges:
        edge.use_freestyle_mark = True

    for mat in obj.material_slots:
        mat.material = bpy.data.materials['diffuse']

    bpy.context.scene.render.filepath = os.path.join(imageDir, name + ".png")
    bpy.ops.render.render(write_still = True)
    bpy.ops.object.delete()
    




        






