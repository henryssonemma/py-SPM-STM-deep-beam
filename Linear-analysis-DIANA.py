import math as math
# Model parameters
folder = "Add-file-path-here"
newProject( folder, 1000, {} )
setModelAnalysisAspects( [ "STRUCT" ] )
setModelDimension( "2D" )
setDefaultMeshOrder( "LINEAR" )
setDefaultMesherType( "HEXQUAD" )
showWorkingPlane( True )
setViewPoint( "FRONT" )
setViewPoint( "BACK" )
setViewPoint( "TOP" )

    # Force
Fy = -3000000

    # Dimensions

# Dimensions concrete beam
cc = 0.05;
L_beam = 4
#L_beam = 3
#L_beam = 3.5
#L_beam = 4.5
#L_beam = 5
#L_beam = 5.5
#L_beam = 6
#L_beam = 6.5
#L_beam = 7
#L_beam = 7.5
#L_beam = 8
#L_beam = 8.5
#L_beam = 9
H_beam = 3
T_beam = 0.400

# Dimensions supports
sup_thic=0.1
sup_length=0.4
left_sup_s=0
right_sup_s=L_beam-left_sup_s-sup_length;

# Dimensions loading plate
plate_length = 0.4;
plate_x = L_beam/2-plate_length/2;
plate_y = H_beam;
plate_thic = 0.1;

# Dimensions hole
#hole_dx = 0.5;
#hole_dx = 0.75;
#hole_dx = 1;
#hole_dx = 1.25;
hole_dx = 1.5;
hole_dy = hole_dx;
hole_x = L_beam/2-hole_dx/2
hole_y = H_beam/2-hole_dy/2

    # Material input
    
# Concrete
# C30/37

# Plates
E_s = 2e+15         # Youngs modulus (fictive value)
v_s = 0.2           # Poissons
p_s = 7800          # Mass density

    # Create materials

#  Fictive plate
addMaterial( "Steel", "MCSTEL", "ISOTRO", [] )
setParameter( "MATERIAL", "Steel", "LINEAR/ELASTI/YOUNG", E_s )
setParameter( "MATERIAL", "Steel", "LINEAR/ELASTI/POISON", v_s )
setParameter( "MATERIAL", "Steel", "LINEAR/MASS/DENSIT", p_s )

# Concrete
addMaterial( "Concrete", "CONCDC", "EN1992", [ "TOTCRK" ] )
setParameter( "MATERIAL", "Concrete", "EC2CON/NORMAL/CLASS", "C30/37" )


    # Create geometries

# Create thickness
addGeometry( "Thickness", "SHEET", "MEMBRA", [] )
setParameter( "GEOMET", "Thickness", "THICK", T_beam )

    # Assign properties and model geometries

# Concrete wall geometry
createSheet( "beam", [ [ 0, 0 ], [ L_beam, 0 ], [ L_beam, H_beam ], [ 0, H_beam ] ] )
setElementClassType( "SHAPE", [ "beam" ], "MEMBRA" )
assignMaterial( "Concrete", "SHAPE", [ "beam" ] )
assignGeometry( "Thickness", "SHAPE", [ "beam" ] )

# Support and load plates
createSheet( "Support_plate_left", [ [ left_sup_s, 0 ], [ left_sup_s+sup_length, 0 ], [ left_sup_s+sup_length, -sup_thic ], [ left_sup_s, -sup_thic ] ] )
createSheet( "Support_plate_right", [ [ right_sup_s, 0 ], [ right_sup_s+sup_length, 0 ], [ right_sup_s+sup_length, -sup_thic ], [ right_sup_s, -sup_thic ] ] )
createSheet( "Load_plate", [ [ plate_x, plate_y ], [ plate_x+plate_length, plate_y ], [ plate_x+plate_length, plate_y+plate_thic ], [ plate_x, plate_y+plate_thic ] ] )
setElementClassType( "SHAPE", [ "Support_plate_left", "Load_plate", "Support_plate_right" ] , "MEMBRA" )
assignMaterial( "Steel", "SHAPE", [ "Support_plate_left", "Load_plate", "Support_plate_right" ]  )
assignGeometry( "Thickness", "SHAPE", [ "Load_plate", "Support_plate_left", "Support_plate_right" ] )

# Create hole
createSheet( "hole", [ [ hole_x, hole_y ], [ hole_x+hole_dx, hole_y ], [hole_x+hole_dx, hole_y+hole_dy ], [ hole_x, hole_y+hole_dy ] ] )
subtract( "beam", [ "hole" ], False, True )

    # Supports, BC

# Create points on support plates
createPointBody( "point 1", [ right_sup_s+sup_length/2, -0.5 ] )
projection( "SHAPEEDGE", "Support_plate_right", [ [ right_sup_s+sup_length/2, -sup_thic ] ], [ "point 1" ], [ 0, 1 ], True )
removeShape( [ "point 1" ] )
createPointBody( "point 1", [ left_sup_s+sup_length/2, -0.5 ] )
projection( "SHAPEEDGE", "Support_plate_left", [ [ left_sup_s+sup_length/2, -sup_thic ] ], [ "point 1" ], [ 0, 1 ], True )
removeShape( [ "point 1" ] )

# Ad BC
addSet( "GEOMETRYSUPPORTSET", "Support set 1" )
createVertexSupport( "Support Left", "Support set 1" )
setParameter( "GEOMETRYSUPPORT", "Support Left", "AXES", [ 1, 2 ] )
setParameter( "GEOMETRYSUPPORT", "Support Left", "TRANSL", [ 1, 1 ] )
setParameter( "GEOMETRYSUPPORT", "Support Left", "ROTATI", [ 0, 0, 0 ] )
attach( "GEOMETRYSUPPORT", "Support Left", "Support_plate_left", [ [ left_sup_s+sup_length/2, -sup_thic ] ] )
createVertexSupport( "Support Right", "Support set 1" )
setParameter( "GEOMETRYSUPPORT", "Support Right", "AXES", [ 1, 2 ] )
setParameter( "GEOMETRYSUPPORT", "Support Right", "TRANSL", [ 0, 1 ] )
setParameter( "GEOMETRYSUPPORT", "Support Right", "ROTATI", [ 0, 0, 0 ] )
attach( "GEOMETRYSUPPORT", "Support Right", "Support_plate_right", [ [ right_sup_s+sup_length/2, -sup_thic ] ] )

    # Load plates, loads

# Create points on load plates
createPointBody( "point 1", [plate_x+plate_length/2, H_beam+0.5 ] )
projection( "SHAPEEDGE", "Load_plate", [ [ plate_x+plate_length/2, plate_y+plate_thic ] ], [ "point 1" ], [ 0, -1 ], True )
removeShape( [ "point 1" ] )

# Add load on load plate
addSet( "GEOMETRYLOADSET", "Load case 1" )
createVertexLoad( "Load 1", "Load case 1" )
setParameter( "GEOMETRYLOAD", "Load 1", "FORCE/VALUE", Fy )
setParameter( "GEOMETRYLOAD", "Load 1", "FORCE/DIRECT", 2 )
attach( "GEOMETRYLOAD", "Load 1", "Load_plate", [ [plate_x+plate_length/2, H_beam+0.5 ] ] )

    # Mesh

# Mesh concrete
setElementSize( [ "beam" ], 0.1, 0.5, True )
setMesherType( [ "beam" ], "HEXQUAD" )

# Mesh plates
setElementSize( [ "Load_plate", "Support_plate_left", "Support_plate_right" ], 0.1, 0.5, True )
clearMesherType( [ "Load_plate", "Support_plate_left", "Support_plate_right" ] )

# Mesh
generateMesh( [] )
hideView( "GEOM" )
showView( "MESH" )

# Linear analysis
addAnalysis( "Liner analysis" )
addAnalysisCommand( "Liner analysis", "LINSTA", "Structural linear static" )

# Run analyses 
runSolver( [ "Liner analysis" ] )
showView( "RESULT" )
selectResult( { "component": "SYY", "result": "Cauchy Total Stresses", "type": "Element", "location": "node" } )
setResultPlot( "tensor-rosette" )
setViewSettingValue( "view setting", "RESULT/TENSOR/SIZE/FACTOR", 10 )
setViewSettingValue( "view setting", "RESULT/TENSOR/SIZE/CRES", 5 )
setViewSettingValue( "view setting", "RESULT/TENSOR/SIZE/MAXMOD", "CLIP" )
setViewSettingValue( "view setting", "RESULT/TENSOR/SIZE/MAXVAL", 20000000 )
setViewSettingValue( "view setting", "RESULT/TENSOR/SIZE/THICK", 2 )
setViewSettingValue( "view setting", "RESULT/TENSOR/EQUIDI/NUMBER", 10 )
setViewSettingValue( "view setting", "RESULT/BCKGRD/OPACIT", 0.5 )
fitAll(  )
