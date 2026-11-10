import pathlib
import bpy

sMaterialName = "PointCloudMaterial"
sAttributeName = "Col"

class ECSVIndex:
    X = 0
    Y = 1
    Z = 2

    R = 3
    G = 4
    B = 5


def SetVertexColor( _obj, _color ):
    vColors =_obj.data.vertex_colors[sAttributeName]

    for idx in range( len( vColors.data)):
        vColors.data[idx].color = _color

    pass


sFilename="C:/Users/Wunstorf/Documents/My Blender/depthField/depth.2022-09-24.csv" 
pathFile = pathlib.Path(sFilename)
sText = pathFile.read_text()
lsText = sText.split( "\n" )

if sMaterialName in bpy.data.materials.keys():
    matPointCloud = bpy.data.materials[sMaterialName]
    bpy.data.materials.remove(matPointCloud)

matPointCloud = bpy.data.materials.new(name=sMaterialName)  # Create a material.
matPointCloud.use_nodes = True
matNodes = matPointCloud.node_tree.nodes
# a new material node tree already has a BSDF and material output node
nodeOutput = matNodes['Material Output']
nodeBSDF = matNodes['Principled BSDF']
nodeAttribute = matNodes.new("ShaderNodeAttribute")
nodeAttribute.attribute_name = sAttributeName
matNodeLinks = matPointCloud.node_tree.links
matNodeLinks.new( nodeBSDF.inputs['Base Color'], nodeAttribute.outputs['Color'] )

for sLine in lsText:
    lsLine = [s.strip() for s in sLine.split( ",")]
    position = ( 
        float(lsLine[ECSVIndex.X]), 
        float(lsLine[ECSVIndex.Y]), 
        float(lsLine[ECSVIndex.Z])
        )  
    color = ( 
        float(lsLine[ECSVIndex.R])/255.0, 
        float(lsLine[ECSVIndex.G])/255.0, 
        float(lsLine[ECSVIndex.B])/255.0,
        1.0
    )
    bpy.ops.mesh.primitive_cube_add( size=0.75, location=position)
    objCube =bpy.context.active_object
    bpy.ops.object.mode_set(mode = 'VERTEX_PAINT')
    objCube.active_material = matPointCloud
    SetVertexColor( objCube, color)

