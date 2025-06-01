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

    # Force and displacement
dy = -0.008

    # Dimensions

# Dimensions concrete
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
hole_dx = 1;
#hole_dx = 1.25;
#hole_dx = 1.5;
hole_dy = hole_dx;
hole_x = L_beam/2-hole_dx/2
hole_y = H_beam/2-hole_dy/2

# Dimensions for tension reinforcements, change area when needed
tie1_area = 0.0018
tie1_d = 0.02
tie1_area_one_bar = math.pi*(tie1_d/2)**2
n_1_bars = tie1_area/tie1_area_one_bar
tie1_perimeter_one = math.pi*tie1_d
tie1_perimeter = tie1_perimeter_one*n_1_bars
tie_1_x1 = 0
tie_1_y1 = 0.08
tie_1_x2 = 4
tie_1_y2 = 0.08

tie2_area = 0.0007
tie2_d = 0.02
tie2_area_one_bar = math.pi*(tie2_d/2)**2
n_2_bars = tie2_area/tie2_area_one_bar
tie2_perimeter_one = math.pi*tie2_d
tie2_perimeter = tie2_perimeter_one*n_2_bars
tie_2_x1 = 1.4
tie_2_y1 = 2.08
tie_2_x2 = 2.6
tie_2_y2 = 2.08

    # Add more ties if needed

    # Grid general
A_min_grid = 400*1e-6

grid_y = 0.2827                                     # y-distance between horizontal bars
grid_x = 0.2827                                     # x-distance between vertical bars
nh_bars = math.floor((H_beam-2*cc)/grid_y);         # number of horizontal bars
nv_bars = math.floor((L_beam-2*cc)/grid_x);         # number of vertical bars
grid_x = round(1000*(L_beam-2*cc)/nv_bars)/1000;    # x-distance between vertical bars
grid_y = round(1000*(H_beam-2*cc)/nh_bars)/1000;    # y-distance between horizontal bars

grid_h_area = 2*A_min_grid*grid_y;              # Area of horizontal bars (dubble because of two sides)
grid_h_d = 2*(((grid_h_area/2)/math.pi)**0.5);  # Diameter of horizontal bars 
grid_h_perimeter = 2*math.pi*grid_h_d;          # Perimeter of horizontal bars (dubble because of two sides)

grid_v_area = 2*A_min_grid*grid_x;              # Area of vertical bars (dubble because of two sides)
grid_v_d = 2*(((grid_v_area/2)/math.pi)**0.5);  # Diameter of horizontal bars
grid_v_perimeter = 2*math.pi*grid_v_d ;         # Perimeter of vertical bars (dubble because of two sides)


    # P4 - panel
A_min_grid_4 = 610*1e-6-400*1e-6
grid_y_4 = 0.2500                                                   # y-distance between horizontal bars
grid_x_4 = 0.2500                                                   # x-distance between vertical bars

x_P4_1 = cc
y_P4_1 = hole_y+grid_y_4/2
x_P4_2 = hole_x-cc
y_P4_2 = hole_y+hole_dy-grid_y_4/2

nh_bars_4 = math.floor((y_P4_2-y_P4_1)/grid_y_4);        # number of horizontal bars
nv_bars_4 = math.floor((x_P4_2-x_P4_1)/grid_x_4);        # number of vertical bars
grid_x_4 = round(1000*(x_P4_2-x_P4_1)/nv_bars_4)/1000;   # x-distance between vertical bars
grid_y_4 = round(1000*(y_P4_2-y_P4_1)/nh_bars_4)/1000;   # y-distance between horizontal bars

grid_4_h_area = 2*A_min_grid_4*grid_y_4;                # Area of horizontal bars (dubble because of two sides)
grid_4_h_d = 2*(((grid_4_h_area/2)/math.pi)**0.5);      # Diameter of horizontal bars 
grid_4_h_perimeter = 2*math.pi*grid_4_h_d;              # Perimeter of horizontal bars (dubble because of two sides)

grid_4_v_area = 2*A_min_grid_4*grid_x_4;                # Area of vertical bars (dubble because of two sides)
grid_4_v_d = 2*(((grid_4_v_area/2)/math.pi)**0.5);      # Diameter of horizontal bars
grid_4_v_perimeter = 2*math.pi*grid_4_v_d ;             # Perimeter of vertical bars (dubble because of two sides)

    # P5 - panel
A_min_grid_5 = 610*1e-6-400*1e-6
grid_y_5 = 0.2500                                                   # y-distance between horizontal bars
grid_x_5 = 0.2500                                                   # x-distance between vertical bars

x_P5_1 = hole_x+hole_dx+cc
y_P5_1 = hole_y+grid_y_5/2
x_P5_2 = L_beam-cc
y_P5_2 = hole_y+hole_dy-grid_y_5/2

nh_bars_5 = math.floor((y_P5_2-y_P5_1)/grid_y_5);        # number of horizontal bars
nv_bars_5 = math.floor((x_P5_2-x_P5_1)/grid_x_5);        # number of vertical bars
grid_x_5 = round(1000*(x_P5_2-x_P5_1)/nv_bars_5)/1000;   # x-distance between vertical bars
grid_y_5 = round(1000*(y_P5_2-y_P5_1)/nh_bars_5)/1000;   # y-distance between horizontal bars

grid_5_h_area = 2*A_min_grid_5*grid_y_5;                # Area of horizontal bars (dubble because of two sides)
grid_5_h_d = 2*(((grid_5_h_area/2)/math.pi)**0.5);      # Diameter of horizontal bars 
grid_5_h_perimeter = 2*math.pi*grid_5_h_d;              # Perimeter of horizontal bars (dubble because of two sides)

grid_5_v_area = 2*A_min_grid_5*grid_x_5;                # Area of vertical bars (dubble because of two sides)
grid_5_v_d = 2*(((grid_5_v_area/2)/math.pi)**0.5);      # Diameter of horizontal bars
grid_5_v_perimeter = 2*math.pi*grid_5_v_d ;             # Perimeter of vertical bars (dubble because of two sides)

    # P7 - panel
A_min_grid_7 = 2034*1e-6-400*1e-6
grid_y_7 = 0.2500                                                   # y-distance between horizontal bars
grid_x_7 = 0.1                                                   # x-distance between vertical bars

x_P7_1 = hole_x+grid_x_7/2
y_P7_1 = hole_y+hole_dy+cc
x_P7_2 = hole_x+hole_dx/2-grid_x_7/2
y_P7_2 = H_beam-cc

nh_bars_7 = math.floor((y_P7_2-y_P7_1)/grid_y_7);        # number of horizontal bars
nv_bars_7 = math.floor((x_P7_2-x_P7_1)/grid_x_7);        # number of vertical bars
grid_x_7 = round(1000*(x_P7_2-x_P7_1)/nv_bars_7)/1000;   # x-distance between vertical bars
grid_y_7 = round(1000*(y_P7_2-y_P7_1)/nh_bars_7)/1000;   # y-distance between horizontal bars

grid_7_h_area = 2*A_min_grid_7*grid_y_7;                # Area of horizontal bars (dubble because of two sides)
grid_7_h_d = 2*(((grid_7_h_area/2)/math.pi)**0.5);      # Diameter of horizontal bars 
grid_7_h_perimeter = 2*math.pi*grid_7_h_d;              # Perimeter of horizontal bars (dubble because of two sides)

grid_7_v_area = 2*A_min_grid_7*grid_x_7;                # Area of vertical bars (dubble because of two sides)
grid_7_v_d = 2*(((grid_7_v_area/2)/math.pi)**0.5);      # Diameter of horizontal bars
grid_7_v_perimeter = 2*math.pi*grid_7_v_d ;             # Perimeter of vertical bars (dubble because of two sides)

    # P8 - panel
A_min_grid_8 = 2034*1e-6-400*1e-6
grid_y_8 = 0.2500                                                   # y-distance between horizontal bars
grid_x_8 = 0.1                                                   # x-distance between vertical bars

x_P8_1 = hole_x+hole_dx/2+grid_x_8/2
y_P8_1 = hole_y+hole_dy+cc
x_P8_2 = hole_x+hole_dx-grid_x_8/2
y_P8_2 = H_beam-cc

nh_bars_8 = math.floor((y_P8_2-y_P8_1)/grid_y_8);        # number of horizontal bars
nv_bars_8 = math.floor((x_P8_2-x_P8_1)/grid_x_8);        # number of vertical bars
grid_x_8 = round(1000*(x_P8_2-x_P8_1)/nv_bars_8)/1000;   # x-distance between vertical bars
grid_y_8 = round(1000*(y_P8_2-y_P8_1)/nh_bars_8)/1000;   # y-distance between horizontal bars

grid_8_h_area = 2*A_min_grid_8*grid_y_8;                # Area of horizontal bars (dubble because of two sides)
grid_8_h_d = 2*(((grid_8_h_area/2)/math.pi)**0.5);      # Diameter of horizontal bars 
grid_8_h_perimeter = 2*math.pi*grid_8_h_d;              # Perimeter of horizontal bars (dubble because of two sides)

grid_8_v_area = 2*A_min_grid_8*grid_x_8;                # Area of vertical bars (dubble because of two sides)
grid_8_v_d = 2*(((grid_8_v_area/2)/math.pi)**0.5);      # Diameter of horizontal bars
grid_8_v_perimeter = 2*math.pi*grid_8_v_d ;             # Perimeter of vertical bars (dubble because of two sides)

    # Add more panels if needed


    # Material input
    
# Concrete
# C30/37
f_ck = 30000000
I_E = 66.96         # Mode I-tensile fracture energy  
e_m = 0.0020        # Strain at maximum stress 
e_u = 0.0105      # Strain at ultimate stress

# Plates
E_s = 2e+11         # Youngs modulus
v_s = 0.2           # Poissons
p_s = 7800          # Mass density


# Reinforcement
E_R = 2e+11             # Youngs modulus

f_yk = 500000000        # Yeild stress
e_ys = f_yk/E_R         # Yeild strain

f_yk1 = 505000000       # Yeild stress 
e_ys1 = 0.02            # Strain plastic

f_u = 600000000         # Ultimate stress
e_us = 0.2              # Ultimate strain

f_u2 = 500000000        # Fail stress
e_us2 = 0.20001         # Fail strain

v_R = 0.2               # Poissons
p_R = 7800              # Mass density
S_N_R = 1e+15           # Normal stiffness modulus (bond-slip)
S_S_R = 1e+12           # Shear Stiffness modulus (bond-slip)
t_m_R = 13.6931e+6      # Maxium shear stress (bond-slip)
t_u_R = 5.4772e+6       # Ultimat shear stress (bond-slip)
s_0 = 0.0001            # Initial slip section (bond-slip)
s_1 = 0.001             # Relative slip section (bond-slip)
s_2 = 0.002             # Relative slip section (bond-slip)
s_3 = 0.005             # Relative slip section (bond-slip)
a_e = 0.4               # Exponent alpa (bond-slip)

# Steel plate, extremly stiff becuse deformation there not interesting
addMaterial( "Steel", "MCSTEL", "ISOTRO", [] )
setParameter( "MATERIAL", "Steel", "LINEAR/ELASTI/YOUNG", E_s )
setParameter( "MATERIAL", "Steel", "LINEAR/ELASTI/POISON", v_s )
setParameter( "MATERIAL", "Steel", "LINEAR/MASS/DENSIT", p_s)

    # Create materials

# Material: concrete
addMaterial( "Concrete", "CONCDC", "EN1992", [ "TOTCRK" ] )
setParameter( "MATERIAL", "Concrete", "EC2CON/NORMAL/CLASS", "C30/37" )
setParameter( "MATERIAL", "Concrete", "TENSIL/TENCRV", "HORDYK" )
setParameter( "MATERIAL", "Concrete", "TENSIL/GF1", I_E )
setParameter( "MATERIAL", "Concrete", "TENSIL/CBSPEC", "ROTS" ) 
setParameter( "MATERIAL", "Concrete", "COMPRS/COMCRV", "EC2" )
setParameter( "MATERIAL", "Concrete", "COMPRS/EPSC1", e_m )
setParameter( "MATERIAL", "Concrete", "COMPRS/EPSCU", e_u )

# Tensile reinforcement and grid reinforcement
addMaterial( "Reinforcement", "REINFO", "REBOND", [] )
setParameter( "MATERIAL", "Reinforcement", "REBARS/ELASTI/YOUNG", E_R )
setParameter( "MATERIAL", "Reinforcement", "REBARS/POISON/POISON", v_R )
setParameter( "MATERIAL", "Reinforcement", "REBARS/MASS/DENSIT", p_R )
setParameter( "MATERIAL", "Reinforcement", "REBARS/PLATYP", "VMISES" )
setParameter( "MATERIAL", "Reinforcement", "REBARS/PLASTI/TRESSH", "EPSSIG" )
setParameter( "MATERIAL", "Reinforcement", "REBARS/PLASTI/EPSSIG", [] )
setParameter( "MATERIAL", "Reinforcement", "REBARS/PLASTI/EPSSIG", [ e_ys, f_yk, e_ys1, f_yk1, e_us, f_u, e_us2, f_u2 ] ) 
setParameter( "MATERIAL", "Reinforcement", "RESLIP/DSNY", S_N_R ) 
setParameter( "MATERIAL", "Reinforcement", "RESLIP/DSSX", S_S_R ) 
setParameter( "MATERIAL", "Reinforcement", "RESLIP/SHFTYP", "BONDS6" )
setParameter( "MATERIAL", "Reinforcement", "RESLIP/BONDS6/SLPVAL", [ t_m_R, t_u_R, s_0, s_1, s_2, s_3, a_e ] )

    # Create geometries

# Create thickness
addGeometry( "Thickness", "SHEET", "MEMBRA", [] )
setParameter( "GEOMET", "Thickness", "THICK", T_beam )

# Create reinforcement, tension
addGeometry( "tie 1", "RELINE", "REBAR", [] )
setParameter( "GEOMET", "tie 1", "REITYP", "REITRU" )
setParameter( "GEOMET", "tie 1", "REITRU/CROSSE", tie1_area )
setParameter( "GEOMET", "tie 1", "REITRU/PERIME", tie1_perimeter )

addGeometry( "tie 2", "RELINE", "REBAR", [] )
setParameter( "GEOMET", "tie 2", "REITYP", "REITRU" )
setParameter( "GEOMET", "tie 2", "REITRU/CROSSE", tie2_area )
setParameter( "GEOMET", "tie 2", "REITRU/PERIME", tie2_perimeter )

#Add more ties if needed

# Create reinforcement, grid
addGeometry( "grid h", "RELINE", "REBAR", [] )
setParameter( "GEOMET", "grid h", "REITYP", "REITRU" )
setParameter( "GEOMET", "grid h", "REITRU/CROSSE", grid_h_area )
setParameter( "GEOMET", "grid h", "REITRU/PERIME", grid_h_perimeter )

addGeometry( "grid v", "RELINE", "REBAR", [] )
setParameter( "GEOMET", "grid v", "REITYP", "REITRU" )
setParameter( "GEOMET", "grid v", "REITRU/CROSSE", grid_v_area )
setParameter( "GEOMET", "grid v", "REITRU/PERIME", grid_v_perimeter )

addGeometry( "grid 4_h", "RELINE", "REBAR", [] )
setParameter( "GEOMET", "grid 4_h", "REITYP", "REITRU" )
setParameter( "GEOMET", "grid 4_h", "REITRU/CROSSE", grid_4_h_area )
setParameter( "GEOMET", "grid 4_h", "REITRU/PERIME", grid_4_h_perimeter )

addGeometry( "grid 4_v", "RELINE", "REBAR", [] )
setParameter( "GEOMET", "grid 4_v", "REITYP", "REITRU" )
setParameter( "GEOMET", "grid 4_v", "REITRU/CROSSE", grid_4_v_area )
setParameter( "GEOMET", "grid 4_v", "REITRU/PERIME", grid_4_v_perimeter )

addGeometry( "grid 5_h", "RELINE", "REBAR", [] )
setParameter( "GEOMET", "grid 5_h", "REITYP", "REITRU" )
setParameter( "GEOMET", "grid 5_h", "REITRU/CROSSE", grid_5_h_area )
setParameter( "GEOMET", "grid 5_h", "REITRU/PERIME", grid_5_h_perimeter )

addGeometry( "grid 5_v", "RELINE", "REBAR", [] )
setParameter( "GEOMET", "grid 5_v", "REITYP", "REITRU" )
setParameter( "GEOMET", "grid 5_v", "REITRU/CROSSE", grid_5_v_area )
setParameter( "GEOMET", "grid 5_v", "REITRU/PERIME", grid_5_v_perimeter )

addGeometry( "grid 7_h", "RELINE", "REBAR", [] )
setParameter( "GEOMET", "grid 7_h", "REITYP", "REITRU" )
setParameter( "GEOMET", "grid 7_h", "REITRU/CROSSE", grid_7_h_area )
setParameter( "GEOMET", "grid 7_h", "REITRU/PERIME", grid_7_h_perimeter )

addGeometry( "grid 7_v", "RELINE", "REBAR", [] )
setParameter( "GEOMET", "grid 7_v", "REITYP", "REITRU" )
setParameter( "GEOMET", "grid 7_v", "REITRU/CROSSE", grid_7_v_area )
setParameter( "GEOMET", "grid 7_v", "REITRU/PERIME", grid_7_v_perimeter )

addGeometry( "grid 8_h", "RELINE", "REBAR", [] )
setParameter( "GEOMET", "grid 8_h", "REITYP", "REITRU" )
setParameter( "GEOMET", "grid 8_h", "REITRU/CROSSE", grid_8_h_area )
setParameter( "GEOMET", "grid 8_h", "REITRU/PERIME", grid_8_h_perimeter )

addGeometry( "grid 8_v", "RELINE", "REBAR", [] )
setParameter( "GEOMET", "grid 8_v", "REITYP", "REITRU" )
setParameter( "GEOMET", "grid 8_v", "REITRU/CROSSE", grid_8_v_area )
setParameter( "GEOMET", "grid 8_v", "REITRU/PERIME", grid_8_v_perimeter )


    # Assign properties and model geometries

# Concrete beam geometry
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

# Reinforcement, tension
addSet( "GEOMETRYREINFOSET", "Tie 1" )
createLine( "Tie_1", [ tie_1_x1, tie_1_y1 ], [ tie_1_x2, tie_1_y2 ] )
setReinforcementType( "GEOMETRYREINFOSET", "Tie 1", "TRUSS_BOND_SLIP" )
assignMaterial( "Reinforcement", "GEOMETRYREINFOSET", [ "Tie 1" ] )
assignGeometry( "tie 1", "GEOMETRYREINFOSET", [ "Tie 1" ] )

addSet( "GEOMETRYREINFOSET", "Tie 2" )
createLine( "Tie_2", [ tie_2_x1, tie_2_y1 ], [ tie_2_x2, tie_2_y2 ] )
setReinforcementType( "GEOMETRYREINFOSET", "Tie 2", "TRUSS_BOND_SLIP" )
assignMaterial( "Reinforcement", "GEOMETRYREINFOSET", [ "Tie 2" ] )
assignGeometry( "tie 2", "GEOMETRYREINFOSET", [ "Tie 2" ] )

# Ad more ties if neccesary

# Reinforcement, grid 

    # general
addSet( "GEOMETRYREINFOSET", "Grid horizontal" )
for i in range(nh_bars+1):
    bar_segment= 'Bar_h' +  str(i) 
    createLine( bar_segment , [ 0 , cc + grid_y*i ], [ L_beam, cc + grid_y*i ] )
    setReinforcementType( "GEOMETRYREINFOSET", "Grid horizontal", "TRUSS_BOND_SLIP" )
    assignMaterial( "Reinforcement", "GEOMETRYREINFOSET", [ "Grid horizontal" ] )
    assignGeometry( "grid h", "GEOMETRYREINFOSET", [ "Grid horizontal" ] )

addSet( "GEOMETRYREINFOSET", "Grid vertical" )
for j in range(nv_bars+1):
    bar_segment= 'Bar_v' +  str(j) 
    createLine( bar_segment , [ cc + grid_x*j, 0 ], [ cc + grid_x*j, H_beam] )
    setReinforcementType( "GEOMETRYREINFOSET", "Grid vertical", "TRUSS_BOND_SLIP" )
    assignMaterial( "Reinforcement", "GEOMETRYREINFOSET", [ "Grid vertical" ] )
    assignGeometry( "grid v", "GEOMETRYREINFOSET", [ "Grid vertical" ] )

 #Take away reinforcement in hole
for i in range(nh_bars_1+1):
    bar_segment= 'Bar_1h' +  str(i) 
    cut( bar_segment, [ "beam" ], True, True )
    removeShape( [bar_segment + '_1' ] )

for i in range(nv_bars_1+1):
    bar_segment= 'Bar_1v' +  str(i) 
    cut( bar_segment, [ "beam" ], True, True )
    removeShape( [bar_segment + '_1' ] )


    #P4
addSet( "GEOMETRYREINFOSET", "Grid 4 horizontal" )
for i in range(nh_bars_4+1):
    bar_segment= 'Bar_4h' +  str(i) 
    createLine( bar_segment , [ x_P4_1 , y_P4_1 + grid_y_4*i ], [ x_P4_2, y_P4_1 + grid_y_4*i ] )
    setReinforcementType( "GEOMETRYREINFOSET", "Grid 4 horizontal", "TRUSS_BOND_SLIP" )
    assignMaterial( "Reinforcement", "GEOMETRYREINFOSET", [ "Grid 4 horizontal" ] )
    assignGeometry( "grid 4_h", "GEOMETRYREINFOSET", [ "Grid 4 horizontal" ] )

addSet( "GEOMETRYREINFOSET", "Grid 4 vertical" )
for j in range(nv_bars_4+1):
    bar_segment= 'Bar_4v' +  str(j) 
    createLine( bar_segment , [ x_P4_1 + grid_x_4*j, y_P4_1 ], [ x_P4_1 + grid_x_4*j, y_P4_2 ] )
    setReinforcementType( "GEOMETRYREINFOSET", "Grid 4 vertical", "TRUSS_BOND_SLIP" )
    assignMaterial( "Reinforcement", "GEOMETRYREINFOSET", [ "Grid 4 vertical" ] )
    assignGeometry( "grid 4_v", "GEOMETRYREINFOSET", [ "Grid 4 vertical" ] )

     #P5
addSet( "GEOMETRYREINFOSET", "Grid 5 horizontal" )
for i in range(nh_bars_5+1):
    bar_segment= 'Bar_5h' +  str(i) 
    createLine( bar_segment , [ x_P5_1 , y_P5_1 + grid_y_5*i ], [ x_P5_2, y_P5_1 + grid_y_5*i ] )
    setReinforcementType( "GEOMETRYREINFOSET", "Grid 5 horizontal", "TRUSS_BOND_SLIP" )
    assignMaterial( "Reinforcement", "GEOMETRYREINFOSET", [ "Grid 5 horizontal" ] )
    assignGeometry( "grid 5_h", "GEOMETRYREINFOSET", [ "Grid 5 horizontal" ] )

addSet( "GEOMETRYREINFOSET", "Grid 5 vertical" )
for j in range(nv_bars_5+1):
    bar_segment= 'Bar_5v' +  str(j) 
    createLine( bar_segment , [ x_P5_1 + grid_x_5*j, y_P5_1 ], [ x_P5_1 + grid_x_5*j, y_P5_2 ] )
    setReinforcementType( "GEOMETRYREINFOSET", "Grid 5 vertical", "TRUSS_BOND_SLIP" )
    assignMaterial( "Reinforcement", "GEOMETRYREINFOSET", [ "Grid 5 vertical" ] )
    assignGeometry( "grid 5_v", "GEOMETRYREINFOSET", [ "Grid 5 vertical" ] )

    #P7
addSet( "GEOMETRYREINFOSET", "Grid 7 horizontal" )
for i in range(nh_bars_7+1):
    bar_segment= 'Bar_7h' +  str(i) 
    createLine( bar_segment , [ x_P7_1 , y_P7_1 + grid_y_7*i ], [ x_P7_2, y_P7_1 + grid_y_7*i ] )
    setReinforcementType( "GEOMETRYREINFOSET", "Grid 7 horizontal", "TRUSS_BOND_SLIP" )
    assignMaterial( "Reinforcement", "GEOMETRYREINFOSET", [ "Grid 7 horizontal" ] )
    assignGeometry( "grid 7_h", "GEOMETRYREINFOSET", [ "Grid 7 horizontal" ] )

addSet( "GEOMETRYREINFOSET", "Grid 7 vertical" )
for j in range(nv_bars_7+1):
    bar_segment= 'Bar_7v' +  str(j) 
    createLine( bar_segment , [ x_P7_1 + grid_x_7*j, y_P7_1 ], [ x_P7_1 + grid_x_7*j, y_P7_2 ] )
    setReinforcementType( "GEOMETRYREINFOSET", "Grid 7 vertical", "TRUSS_BOND_SLIP" )
    assignMaterial( "Reinforcement", "GEOMETRYREINFOSET", [ "Grid 7 vertical" ] )
    assignGeometry( "grid 7_v", "GEOMETRYREINFOSET", [ "Grid 7 vertical" ] )

    #P8
addSet( "GEOMETRYREINFOSET", "Grid 8 horizontal" )
for i in range(nh_bars_8+1):
    bar_segment= 'Bar_8h' +  str(i) 
    createLine( bar_segment , [ x_P8_1 , y_P8_1 + grid_y_8*i ], [ x_P8_2, y_P8_1 + grid_y_8*i ] )
    setReinforcementType( "GEOMETRYREINFOSET", "Grid 8 horizontal", "TRUSS_BOND_SLIP" )
    assignMaterial( "Reinforcement", "GEOMETRYREINFOSET", [ "Grid 8 horizontal" ] )
    assignGeometry( "grid 8_h", "GEOMETRYREINFOSET", [ "Grid 8 horizontal" ] )

addSet( "GEOMETRYREINFOSET", "Grid 8 vertical" )
for j in range(nv_bars_8+1):
    bar_segment= 'Bar_8v' +  str(j) 
    createLine( bar_segment , [ x_P8_1 + grid_x_8*j, y_P8_1 ], [ x_P8_1 + grid_x_8*j, y_P8_2 ] )
    setReinforcementType( "GEOMETRYREINFOSET", "Grid 8 vertical", "TRUSS_BOND_SLIP" )
    assignMaterial( "Reinforcement", "GEOMETRYREINFOSET", [ "Grid 8 vertical" ] )
    assignGeometry( "grid 8_v", "GEOMETRYREINFOSET", [ "Grid 8 vertical" ] )

    # Supports and loads, BC, 

# Create points on support plates
createPointBody( "point 1", [ right_sup_s+sup_length/2, -sup_thic ] )
projection( "SHAPEEDGE", "Support_plate_right", [ [ right_sup_s+sup_length/2, -sup_thic ] ], [ "point 1" ], [ 0, 1 ], True )
removeShape( [ "point 1" ] )
createPointBody( "point 1", [ left_sup_s+sup_length/2, -sup_thic ] )
projection( "SHAPEEDGE", "Support_plate_left", [ [ left_sup_s+sup_length/2, -sup_thic ] ], [ "point 1" ], [ 0, 1 ], True )
removeShape( [ "point 1" ] )

# Create points on load plates
createPointBody( "point 1", [plate_x+plate_length/2, H_beam+plate_thic ] )
projection( "SHAPEEDGE", "Load_plate", [ [ plate_x+plate_length/2, plate_y+plate_thic ] ], [ "point 1" ], [ 0, -1 ], True )
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

# Add displacement on load plate
createVertexSupport( "Load plate", "Support set 1" )
setParameter( "GEOMETRYSUPPORT", "Load plate", "AXES", [ 1, 2 ] )
setParameter( "GEOMETRYSUPPORT", "Load plate", "TRANSL", [ 0, 1 ] )
setParameter( "GEOMETRYSUPPORT", "Load plate", "ROTATI", [ 0, 0, 0 ] )
attach( "GEOMETRYSUPPORT", "Load plate", "Load_plate", [ [plate_x+plate_length/2, H_beam+plate_thic ] ] )
addSet( "GEOMETRYLOADSET", "Load case 1" )
createVertexLoad( "Displacement", "Load case 1" )
setParameter( "GEOMETRYLOAD", "Displacement", "LODTYP", "DEFORM" )
setParameter( "GEOMETRYLOAD", "Displacement", "DEFORM/SUPP", "Support 1" )
setParameter( "GEOMETRYLOAD", "Displacement", "DEFORM/TR/VALUE", dy )
setParameter( "GEOMETRYLOAD", "Displacement", "DEFORM/TR/DIRECT", 2 )
attachTo( "GEOMETRYLOAD", "Displacement", "DEFORM/TARGET", "Load_plate", [ [plate_x+plate_length/2, H_beam+plate_thic ] ] )

# Tyings
addSet( "GEOMETRYTYINGSET", "Tying load plate" )
createEdgeTying( "Tying load plate", "Tying load plate" )
setParameter( "GEOMETRYTYING", "Tying load plate", "AXES", [ 1, 2 ] )
setParameter( "GEOMETRYTYING", "Tying load plate", "TRANSL", [ 1, 1 ] )
setParameter( "GEOMETRYTYING", "Tying load plate", "ROTATI", [ 0, 0, 1 ] )
attachTo( "GEOMETRYTYING", "Tying load plate", "MASTER", "Load_plate", [ [plate_x+plate_length/2, H_beam+plate_thic ] ] )

    # Mesh
    
# Geometries for different mesh
createSheet( "loadplate", [ [ plate_x-0.2, plate_y ], [ plate_x+plate_length+0.2, plate_y ], [ plate_x+plate_length+0.2, plate_y-0.1 ], [ plate_x-0.2, plate_y-0.1 ] ] )
createSheet( "supp_rigth", [ [ right_sup_s-0.2, 0 ], [ right_sup_s+sup_length, 0 ], [ right_sup_s+sup_length, 0+0.1 ], [ right_sup_s-0.2, 0+0.1 ] ] )
createSheet( "supp_left", [ [ left_sup_s, 0 ], [ left_sup_s+sup_length+0.2, 0 ], [ left_sup_s+sup_length+0.2, 0+0.1 ], [ left_sup_s, 0+0.1 ] ] )
cut( "beam", [ "loadplate", "supp_rigth", "supp_left" ], False, True )

# Mesh plates
setElementSize( [ "Load_plate", "Support_plate_left", "Support_plate_right" ], 0.1, 0.5, True )
clearMesherType( [ "Load_plate", "Support_plate_left", "Support_plate_right" ] )

# Mesh concrete, large mesh - change if needed
setElementSize( [ "beam_1", "beam_2", "beam_3" ], 0.1, 0.5, True )
clearMesherType( [ "beam_1", "beam_2", "beam_3" ] )

# Mesh concrete
setElementSize( [ "beam" ], 0.1, 0.5, True )
setMesherType( [ "beam" ], "HEXQUAD" )

# Mesh
generateMesh( [] )
hideView( "GEOM" )
showView( "MESH" )

# Linear analysis
addAnalysis( "Liner analysis" )
addAnalysisCommand( "Liner analysis", "LINSTA", "Structural linear static" )

# Non - linear analysis
addAnalysis( "Non-linear analysis" )
addAnalysisCommand( "Non-linear analysis", "NONLIN", "Structural nonlinear" )
setAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "EXECUT(1)/LOAD/STEPS/STEPTY", "EXPLIC" )
setAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "EXECUT(1)/LOAD/STEPS/EXPLIC/SIZES", "0.01000(100)" )
setAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "EXECUT(1)/LOAD/STEPS/AUTOMA/CUTBCK", 0.25 )
setAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "EXECUT(1)/LOAD/STEPS/AUTOMA/SIZES", 5 )
setAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "EXECUT(1)/LOAD/STEPS/AUTOMA/SIZES", 100 )
setAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "EXECUT(1)/ITERAT/MAXITE", 500 )
addAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "EXECUT(1)/ITERAT/CONVER/ENERGY" )
setAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "EXECUT(1)/ITERAT/CONVER/ENERGY", True )

    # All pre added outputs
addAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "OUTPUT(1)" )
setAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "OUTPUT(1)/FILE", "all_primaries" )

    # Cracks
setAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "OUTPUT(5)/SELTYP", "USER" )
setAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "OUTPUT(5)/FILE", "Cracks" )
addAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "OUTPUT(5)/USER/STRESS(1)/CRACK/CAUCHY/LOCAL" )
addAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "OUTPUT(5)/USER/STRAIN(1)/CRACK/GREEN" )
addAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "OUTPUT(5)/USER/STRAIN(2)/CRKSUM/GREEN/LOCAL" )
addAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "OUTPUT(5)/USER/STRAIN(3)/CRKSUM/GREEN/GLOBAL" )
addAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "OUTPUT(5)/USER/STRAIN(4)/CRKSUM/GREEN/PRINCI" )
addAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "OUTPUT(5)/USER/STRAIN(5)/CRKSUM/GREEN/VONMIS" )
addAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "OUTPUT(5)/USER/STRAIN(6)/CRKWDT/GREEN/GLOBAL" )
addAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "OUTPUT(5)/USER/STRAIN(7)/CRKWDT/GREEN/LOCAL" )
addAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "OUTPUT(5)/USER/STRAIN(8)/CRKWDT/GREEN/VONMIS" )
addAnalysisCommandDetail( "Non-linear analysis", "Structural nonlinear", "OUTPUT(5)/USER/STRAIN(9)/CRKWDT/GREEN/PRINCI" )

# Run analyses 
#runSolver( [ "Liner analysis" ] )
#runSolver( [ "Non-linear analysis" ] )

