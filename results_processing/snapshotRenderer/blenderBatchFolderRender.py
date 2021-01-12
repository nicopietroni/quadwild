import bpy
import glob, os, sys

# takes directory basename and returns names of files to be rendered
def getToRenderNames(dirName):
    return [ dirName + ".obj", dirName + "_rem_p0_0_quadrangulation_smooth.obj", dirName + "_rem_p0_0_quadrangulation.obj"]

argv = sys.argv
argv = argv[argv.index("--") + 1:]

baseDirectory = os.path.abspath(argv[0])

print(baseDirectory)

scene = bpy.context.scene
data  = bpy.data

# get all the subdirectories names
subDirectories = next(os.walk(baseDirectory))[1]

for dirName in subDirectories:
    
    # get path to the directory
    dirPath = os.path.join(baseDirectory, dirName)

    # move to subdirectory
    os.chdir(dirPath)    

    meshes = getToRenderNames(dirName)

    for mesh in meshes:
        
        success = bpy.ops.import_scene.obj(filepath=os.path.join('.', mesh), axis_forward='-Z', axis_up='Y')

        # the imported objects...
        obj_objects = bpy.context.selected_objects[:]

        if (len(obj_objects) == 0):
            continue

        print(mesh)

        obj = bpy.context.selected_objects[0]

        bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS', center='MEDIAN')
        obj.location = (0, 0, 0)
        obj.scale = obj.scale / max(obj.dimensions)

        for edge in obj.data.edges:
            edge.use_freestyle_mark = True

        for mat in obj.material_slots:
            mat.material = bpy.data.materials['diffuse']

        bpy.context.scene.render.image_settings.file_format='JPEG'
        bpy.context.scene.render.image_settings.quality=85
        bpy.context.scene.render.filepath = os.path.join(os.getcwd(), mesh + ".jpg")
        bpy.ops.render.render(write_still = True)
        bpy.ops.object.delete()



    









