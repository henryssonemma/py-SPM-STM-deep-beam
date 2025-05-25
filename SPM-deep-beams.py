# Python script for designing and analysing concrete deep beams with the Stringer-Panel Method.
# This script was written for the purpose of the master's thesis titled "Assessing the Stringer-Panel Method as an Alternative Design Approach for Deep Beams - A Comparative Study with the Strut-and-Tie Method", by Alice Hedenberg and Emma Henrysson at the Department of Architecture and Civil Engineering at Chalmers University of Technology during the spring of 2025.
# This script is for designing and analytically assess concrete deep beams
# This code is not universally applicable on all geometries of deep beams. The script has only been designed for and verified for the deep beam geometries andn load cases used in the thesis.
# Parts of the script is borrowed from the script SPM.py written by P.C.J. Hoogenboom, those parts will be clearly marked within this script. Those parts are also slightly modified in order to fit with the format of the geometry vectors created by the first part of the script, but nu changes to the calculations were made.

import numpy as np
from math import floor

# Input data

# Geometry
width = 4 # Width of beam [m]
height = 3 # Height of beam  [m]
c_edge = 0.08 # Distance between concrete edge and center of stringer [m] 
c_supports = 0.2 # Distance between concrete edge and center of stringer aligned with support [m] 
t = 0.4 # Thickness of beam [m]

# Opening geometry 
hole = "yes" # yes or no 
width_hole = 1 # Width of opening [m]
height_hole = 1 # Height of opening [m]
x_bottomleft_hole = (width-width_hole)/2 # Distance in x-direction from the left concrete edge to the bottom left corner of the opening [m]
y_bottomleft_hole = (height-height_hole)/2 # Distance in y-direction from the bottom concrete edge to the bottom left corner of the opening [m]

# Material properties 
nu = 0.2 # Possion's ratio of concrete
E = 32.8e9 # Modulus of elasticity of concrete [Pa]
G = E/(2*(1+nu)) # Shear modulus for concrete [Pa]
fyk = 500000000; # Characteristic yield strength of reinforcement [Pa]

# External forces 
NrForces = 1 # Number of external forces
Force_dim = -3000000 # Design load [N]
Force = -3000000 # Applied load [N] (If a check at any other load than the design load is desired) 
x_force = width/2 # X-coordinate for the external load [m]

# Fixed displacements
NrFixedDofs = 3 # Number of fixed displacements 
Disp = [0,0,0] # Fixed displacements [m]


# Creating stringer- and panel-geometries, degrees of freedom and placements

# Stringers

# Number of stringers
if hole == "yes":
    Str_vert_hole = 2 # Extra stringers because of opening 
    Str_horiz_hole = 2
else:
    Str_vert_hole = 0
    Str_horiz_hole = 0

NrStrs_vert_full= 2 + NrForces + Str_vert_hole # Number of edge-edge vertical stringers 
NrStrs_horiz_full = 2 + Str_horiz_hole # Number of edge-edge horizontal stringers 

# Placement of stringers
if hole == "yes":
    x_left_hole = x_bottomleft_hole-c_edge
    x_right_hole = x_bottomleft_hole+c_edge+width_hole
    x_Strs = np.array([c_supports,x_force,width-c_supports]) 
    x_Strs = np.append(x_Strs, [x_left_hole, x_right_hole])
    x_Strs_vert = np.sort(x_Strs) # Placement in x-direction for vertical stringers 
  
    y_top_hole = y_bottomleft_hole+c_edge+height_hole
    y_bottom_hole = y_bottomleft_hole-c_edge
    y_Strs = np.array([c_edge,height-c_edge]) 
    y_Strs = np.append(y_Strs, [y_top_hole, y_bottom_hole])
    y_Strs_horiz = np.sort(y_Strs) # Placement in y-direction for horizontal stringers 
else:
    x_Strs_vert = np.array([c_supports,x_force,width-c_supports]) # Placement in x-direction for vertical stringers 
    y_Strs_horiz = np.array([c_edge,height-c_edge]) # Placement in y-direction for horizontal stringers 

# Number of stringer segments 
if hole =="yes":
    NrStrs_vert_inside_hole = sum(
    x_left_hole < x < x_right_hole for x in x_Strs_vert
)
    GlobalStrNo_vert_inside_hole = [i for i, x in enumerate(x_Strs_vert) if x_left_hole < x < x_right_hole]

    if NrStrs_vert_inside_hole >= 1:
        NrStr_vert_segments_tot = NrStrs_vert_full*(NrStrs_horiz_full-1)-NrStrs_vert_inside_hole*2 # Total number of vertical stringer segments (assumes that two stringer segments always dissapear per stringer that falls within the opening) 
        NrStr_vert_segments = int((NrStr_vert_segments_tot+NrStrs_vert_inside_hole*2)/NrStrs_vert_full) # Number of segments each vertical edge-edge stringer is divided into, including those where there is no stringer inside and under the opening 
    else:
        NrStr_vert_segments_tot = NrStrs_vert_full*(NrStrs_horiz_full-1) # Total number of vertical stringer segments
        NrStr_vert_segments = int(NrStr_vert_segments_tot/NrStrs_vert_full) # Number of segments each vertical edge-edge stringer is divided into

    NrStr_horiz_segments_tot = NrStrs_horiz_full*(NrStrs_vert_full-1)-2*NrStrs_vert_inside_hole # Total number of horizontal stringer segments, assuming there are always two horizontal stringers below the window, where two stringer segments are "combined" 
    NrStr_horiz_segments = int((NrStr_horiz_segments_tot+2*NrStrs_vert_inside_hole)/NrStrs_horiz_full) # Number of segments each horizontal edge-edge stringer is divided into
else:
    GlobalStrNo_vert_inside_hole = np.nan

    NrStrs_vert_inside_hole = 0

    NrStr_vert_segments_tot = NrStrs_vert_full*(NrStrs_horiz_full-1) # Total number of vertical stringer segments 
    NrStr_vert_segments = int(NrStr_vert_segments_tot/NrStrs_vert_full) # Number of segments each vertical edge-edge stringer is divided into

    NrStr_horiz_segments_tot = NrStrs_horiz_full*(NrStrs_vert_full-1) # Total number of horizontal stringer segments 
    NrStr_horiz_segments = int(NrStr_horiz_segments_tot/NrStrs_horiz_full) # Number of segments each horizontal edge-edge stringer is divided into

NrStrs = NrStr_horiz_segments_tot+NrStr_vert_segments_tot # Total number of stringer segments 

# Lengths of stringer segments 
LengthStr_vert = np.zeros((NrStrs_vert_full, NrStr_vert_segments)) 
for i in range(NrStrs_vert_full):
    for j in range(NrStr_vert_segments):
        length_v = y_Strs_horiz[j+1]-y_Strs_horiz[j]
        LengthStr_vert[i,j] += length_v # (Each row represents one endge-edge stringer)

if hole == "yes":
    for i in range(NrStrs_vert_inside_hole):
        j = GlobalStrNo_vert_inside_hole[i]
        LengthStr_vert[j, :-1] = 0 # Assumes that only one stringer segment is kept above the opening 

LengthStr_vert_vector = LengthStr_vert.flatten() 
LengthStr_vert_vector = LengthStr_vert_vector[LengthStr_vert_vector != 0] 

LengthStr_horiz = np.zeros((NrStrs_horiz_full, NrStr_horiz_segments))
for i in range(NrStrs_horiz_full):
    for j in range(NrStr_horiz_segments):
        length_h = x_Strs_vert[j+1]-x_Strs_vert[j]

        if hole == "yes" and i < 2:
            for vert_str_idx in GlobalStrNo_vert_inside_hole:
                if j == vert_str_idx - 1:  
                    LengthStr_horiz[i, j] += length_h  
                elif j == vert_str_idx:  
                    LengthStr_horiz[i, j - 1] += length_h  
                    LengthStr_horiz[i, j] = 0 
                else:
                    LengthStr_horiz[i, j] += length_h  
        else:
            LengthStr_horiz[i, j] += length_h  

LengthStr_horiz_vector = LengthStr_horiz.flatten()
LengthStr_horiz_vector = LengthStr_horiz_vector[LengthStr_horiz_vector != 0]

Strl = np.concatenate((LengthStr_horiz_vector, LengthStr_vert_vector)) # Lengths of all stringer segments

# Nodal coordinates
NrNodes = NrStrs_horiz_full*NrStrs_vert_full 

x_nodes_strs = np.zeros([NrNodes])
y_nodes_strs = np.zeros([NrNodes])

for i in range(NrStrs_horiz_full):
    for j in range(NrStrs_vert_full):
        index = i * NrStrs_vert_full + j
        x_nodes_strs[index] = x_Strs_vert[j]
        y_nodes_strs[index] = y_Strs_horiz[i]

if hole == "yes": # Remove the nodes for the nodes on vertical stringers that falls inside the opening
    for vert_str_idx in GlobalStrNo_vert_inside_hole:
        for i in range(2): 
            index = (i * NrStrs_vert_full) + vert_str_idx  
            x_nodes_strs[index] = np.nan  
            y_nodes_strs[index] = np.nan  
x_nodes_strs = x_nodes_strs[~np.isnan(x_nodes_strs)]
y_nodes_strs = y_nodes_strs[~np.isnan(y_nodes_strs)]

NrNodes = len(x_nodes_strs) # Total number of nodes

# Degrees of freedom for stringers 
NrStr_horiz_segments_on_row = np.full(NrStrs_horiz_full, NrStr_horiz_segments, dtype=int) # Number of stringer segments for each horizontal stringer 
NrStr_vert_segments_on_row = np.full(NrStrs_vert_full, NrStr_vert_segments, dtype=int) # Number of stringer segments for each vertical stringer 
 
if hole == "yes": 
    NrStr_horiz_segments_on_row[:2] -= 1 # Assumes that there is one stringes segment less on for the two first horizontal stringers
    NrStr_vert_segments_on_row[GlobalStrNo_vert_inside_hole] = 1 # Assumes that there is always one stringer segment kept above the opening for vertical stringer that falls inside the opening

NrDofs = NrNodes*2+NrStr_vert_segments_tot+NrStr_horiz_segments_tot # Number of DOFs

x_dof_matrix = np.zeros((NrStr_horiz_segments_tot, 3), dtype=int) # Degrees of freedom in the x-direction 
dof_counter = 0  
for i in range(NrStrs_horiz_full): 
    for j in range(NrStr_horiz_segments_on_row[i]):
        row_index = np.sum(NrStr_horiz_segments_on_row[:i]) + j 
        if j == 0 and i > 0:
            dof_counter += 1  
            x_dof_matrix[row_index, 0] = dof_counter 
            dof_counter += 1  
        else:
            x_dof_matrix[row_index, 0] = dof_counter
            dof_counter += 1  
        x_dof_matrix[row_index, 1] = dof_counter
        dof_counter += 1  
        x_dof_matrix[row_index, 2] = dof_counter

y_dof_matrix = np.zeros((NrStr_vert_segments_tot, 3), dtype=int) # Degrees of freedom in the y-direction 
dof_counter += 1  
for i in range(NrStrs_vert_full):
    for j in range(NrStr_vert_segments_on_row[i]):
        row_index = np.sum(NrStr_vert_segments_on_row[:i]) + j 
        if j == 0 and i > 0:
            dof_counter += 1  
            y_dof_matrix[row_index, 0] = dof_counter
            dof_counter += 1  
        else:
            y_dof_matrix[row_index, 0] = dof_counter
            dof_counter += 1 
        y_dof_matrix[row_index, 1] = dof_counter
        dof_counter += 1  
        y_dof_matrix[row_index, 2] = dof_counter

StrDof = np.vstack((x_dof_matrix, y_dof_matrix)) # Degrees of freedom for stringer segments


# Panels

NrPanels = NrStr_horiz_segments*NrStr_vert_segments # Number of panels 
if hole == "yes": # Remove panels for the opening, and udner the opening where panels are "combined" 
    NrPanels -= (NrStrs_vert_inside_hole+1 + NrStrs_vert_inside_hole)  

# Degrees of freedom for panels 
panel_dof_matrix = np.zeros((NrPanels, 4), dtype=int)
panel_index = 0 
if hole == "yes":
    cumsum_horiz_segments = np.insert(np.cumsum(NrStr_horiz_segments_on_row), 0, 0)
    cumsum_vert_segments = np.insert(np.cumsum(NrStr_vert_segments_on_row), 0, 0)
    for i in range(NrStrs_horiz_full - 1):  
        NrPanels_on_row = NrStr_horiz_segments_on_row[i]  
        valid_vert_str_indices = [
            j for j in range(NrStrs_vert_full) if j not in GlobalStrNo_vert_inside_hole or i == NrStrs_horiz_full - 2
        ]
        for j_idx, j in enumerate(valid_vert_str_indices[:-1]):  
            if hole == "yes" and i == 1 and j+1 in GlobalStrNo_vert_inside_hole:  
                continue
            if i == NrStrs_horiz_full - 2:  
                j_right = j + 1  
            elif j+1 in GlobalStrNo_vert_inside_hole:
                j_right = j + 2  
            else:
                j_right = j + 1  

            if i == NrStrs_horiz_full - 2 and j+1 in GlobalStrNo_vert_inside_hole: 
                right_dof = y_dof_matrix[cumsum_vert_segments[j_right], 1]  # Right DOF
            else:
                right_dof = y_dof_matrix[cumsum_vert_segments[j_right] + i, 1]  # Right DOF

            if i == NrStrs_horiz_full - 2 and j in GlobalStrNo_vert_inside_hole: 
                left_dof = y_dof_matrix[cumsum_vert_segments[j], 1]  # Left DOF
            else:
                left_dof = y_dof_matrix[cumsum_vert_segments[j] + i, 1]  # Left DOF

            bottom_dof = x_dof_matrix[cumsum_horiz_segments[i] + j_idx, 1]  # Bottom DOF

            if i == 1 and j-1 in GlobalStrNo_vert_inside_hole: 
                top_dof = x_dof_matrix[cumsum_horiz_segments[i + 1] + (j_idx + 1), 1]  # Top DOF
            else:
                top_dof = x_dof_matrix[cumsum_horiz_segments[i + 1] + j_idx, 1]  # Top DOF

            panel_dof_matrix[panel_index] = [bottom_dof, top_dof, left_dof, right_dof] # DOFs for each panel
            panel_index += 1
else: 
    for i in range(NrStrs_horiz_full - 1):  
        for j in range(NrStrs_vert_full - 1):  
            left_dof = y_dof_matrix[j * NrStr_vert_segments + i, 1]  # Left DOF
            right_dof = y_dof_matrix[j * NrStr_vert_segments + i + 1, 1]  # Right DOF
            bottom_dof = x_dof_matrix[i * NrStr_horiz_segments + j, 1]  # Bottom DOF
            top_dof = x_dof_matrix[i * NrStr_horiz_segments + j + NrStr_horiz_segments, 1]  # Top DOF
            panel_dof_matrix[panel_index] = [bottom_dof, top_dof, left_dof, right_dof] # DOFs for each panel
            panel_index += 1

PanelDof = panel_dof_matrix # Panel DOFs for each panel

# Sizes of panels
width_panels = np.zeros(NrPanels)
height_panels = np.zeros(NrPanels)
panel_index = 0  # Reset panel counter
_, nodes_on_row = np.unique(y_nodes_strs, return_counts=True)
cumsum_nodes_rows = np.insert(np.cumsum(nodes_on_row),0,0)

for i in range(NrStrs_horiz_full-1):
    for j in range(nodes_on_row[i]-1):
        if nodes_on_row[i] < nodes_on_row[i+1] and j+1 in GlobalStrNo_vert_inside_hole:
            continue
        elif nodes_on_row[i] < nodes_on_row[i+1] and j+1 > GlobalStrNo_vert_inside_hole[0]:
            j_top_y = j+2
            topleft_y = y_nodes_strs[cumsum_nodes_rows[i+1]+j_top_y]
        else:
            topleft_y = y_nodes_strs[cumsum_nodes_rows[i+1]+j]
        
        bottomleft_y = y_nodes_strs[cumsum_nodes_rows[i]+j]
        bottomleft_x = x_nodes_strs[cumsum_nodes_rows[i]+j]
        bottomright_x = x_nodes_strs[cumsum_nodes_rows[i]+j+1]

        w_panel = bottomright_x-bottomleft_x
        h_panel = topleft_y-bottomleft_y

        width_panels[panel_index] = w_panel
        height_panels[panel_index] = h_panel
        panel_index += 1

Panela = width_panels # Widths of panels
Panelb = height_panels # Heights of panels


# Stiffnesses

# Stringer stiffness
Str_idx = 0
b_store = np.zeros(NrStrs)
StrEA = np.zeros(NrStrs)
_, nodes_on_vertical_row = np.unique(x_nodes_strs, return_counts=True)
cumsum_nodes_vertical_rows = np.insert(np.cumsum(nodes_on_vertical_row),0,0)

for i in range(NrStrs_horiz_full):
    for j in range(nodes_on_row[i]-1):
        if i == 0:
            b = (y_nodes_strs[cumsum_nodes_rows[i+1]+j]-y_nodes_strs[cumsum_nodes_rows[i]+j])/2+y_nodes_strs[cumsum_nodes_rows[i]+j]
        elif i+1 == NrStrs_horiz_full:
            b = height-y_nodes_strs[cumsum_nodes_rows[i]+j]+(y_nodes_strs[cumsum_nodes_rows[i]+j]-y_nodes_strs[cumsum_nodes_rows[i-1]+j])/2
        elif hole == "yes" and j+1 == GlobalStrNo_vert_inside_hole[0] and nodes_on_row[i]<nodes_on_row[i+1]:
            b = c_edge + (y_nodes_strs[cumsum_nodes_rows[i]+j]-y_nodes_strs[cumsum_nodes_rows[i-1]+j])/2
        elif hole == "yes" and (nodes_on_row[i]>nodes_on_row[i-1]) and ((j+1 == GlobalStrNo_vert_inside_hole[0]) or (j == GlobalStrNo_vert_inside_hole[0])):
            b = c_edge + (y_nodes_strs[cumsum_nodes_rows[i+1]+j]-y_nodes_strs[cumsum_nodes_rows[i]+j])/2
        else:
            b = (y_nodes_strs[cumsum_nodes_rows[i]+j]-y_nodes_strs[cumsum_nodes_rows[i-1]+j])/2+(y_nodes_strs[cumsum_nodes_rows[i+1]+j]-y_nodes_strs[cumsum_nodes_rows[i]+j])/2
        
        b_store[Str_idx] = b
        EA = E*t*b
        StrEA[Str_idx] = EA # Extorsional stiffness for each stringer segment
        Str_idx += 1

for i in range(NrStrs_vert_full):
    for j in range(nodes_on_vertical_row[i]-1):
        if i == 2:
            j += nodes_on_vertical_row[0]-2

        if i == 0:
            b = x_nodes_strs[cumsum_nodes_rows[j]+i] + (x_nodes_strs[cumsum_nodes_rows[j]+i+1]-x_nodes_strs[cumsum_nodes_rows[j]+i])/2
        elif i == NrStrs_vert_full-1:
            b = width - x_nodes_strs[cumsum_nodes_rows[j]+NrStr_horiz_segments_on_row[j]] + (x_nodes_strs[cumsum_nodes_rows[j]+NrStr_horiz_segments_on_row[j]]-x_nodes_strs[cumsum_nodes_rows[j]+NrStr_horiz_segments_on_row[j]-1])/2
        elif hole == "yes" and i+1 == GlobalStrNo_vert_inside_hole[0] and j == 1:
            b = c_edge + (x_nodes_strs[cumsum_nodes_rows[j]+i]-x_nodes_strs[cumsum_nodes_rows[j]+i-1])/2
        elif hole == "yes" and i-1 == GlobalStrNo_vert_inside_hole[0] and j == 1:
            b = c_edge + (x_nodes_strs[cumsum_nodes_rows[j]+GlobalStrNo_vert_inside_hole[0]+1]-x_nodes_strs[cumsum_nodes_rows[j]+GlobalStrNo_vert_inside_hole[0]])/2
        elif hole == "yes" and i+1 == GlobalStrNo_vert_inside_hole[0] and j == 2:
            b = (x_nodes_strs[cumsum_nodes_rows[j]+i]-x_nodes_strs[cumsum_nodes_rows[j]+i-1])/2 + (x_nodes_strs[cumsum_nodes_rows[j]+i+1]-x_nodes_strs[cumsum_nodes_rows[j]+i])/2
        elif hole == "yes" and i-1 == GlobalStrNo_vert_inside_hole[0] and j == 2:
            b = (x_nodes_strs[cumsum_nodes_rows[j]+i+1]-x_nodes_strs[cumsum_nodes_rows[j]+i])/2 + (x_nodes_strs[cumsum_nodes_rows[j]+i]-x_nodes_strs[cumsum_nodes_rows[j]+i-1])/2
        elif hole == "yes" and i-1 == GlobalStrNo_vert_inside_hole[0] and j == 0:
            i = 2
            b = (x_nodes_strs[cumsum_nodes_rows[j]+i]-x_nodes_strs[cumsum_nodes_rows[j]+i-1])/2 + (x_nodes_strs[cumsum_nodes_rows[j]+i+1]-x_nodes_strs[cumsum_nodes_rows[j]+i])/2
            i = 3
        else:
            b = (x_nodes_strs[cumsum_nodes_rows[j]+i]-x_nodes_strs[cumsum_nodes_rows[j]+i-1])/2 + (x_nodes_strs[cumsum_nodes_rows[j]+i+1]-x_nodes_strs[cumsum_nodes_rows[j]+i])/2

        b_store[Str_idx] = b
        EA = E*t*b
        StrEA[Str_idx] = EA # Extorsional stiffness for each stringer segment
        Str_idx += 1

# Panel shear stiffness
PanelGt = np.full(NrPanels,G*t)

# Fixed dofs
FixedDof = np.zeros(NrFixedDofs,dtype=int)
for i in range(NrFixedDofs):
    if i == 0:
        FixedDof[i] = x_dof_matrix[0,0]
    if i == 1:
        FixedDof[i] = y_dof_matrix[0,0]
    if i == 2:
        FixedDof[i] = y_dof_matrix[NrStr_vert_segments_tot-NrStr_vert_segments,0]


# Calculations for reinforcement
# -- The following calculations are taken from the script SPM.py by Hoogenboom (with some additions) --

# Initialize stiffness matrix
s = np.zeros((NrDofs, NrDofs))

# Assemble stringers
sl = np.array([[4, -6, 2], [-6, 12, -6], [2, -6, 4]])
for i in range(NrStrs):
    c = StrEA[i] / Strl[i]
    for k in range(3):
        kk = StrDof[i, k]
        for l in range(3):
            ll = StrDof[i, l]
            s[kk, ll] += c * sl[k, l]

# Assemble panels
for i in range(NrPanels):
    a = Panela[i] / Panelb[i]
    e = 1 / a
    c = PanelGt[i]
    sl = np.array([[a, -a, 1, -1], [-a, a, -1, 1], [1, -1, e, -e], [-1, 1, -e, e]])
    for k in range(4):
        kk = PanelDof[i, k] 
        for l in range(4):
            ll = PanelDof[i, l] 
            s[kk, ll] += c * sl[k, l]

# Process imposed forces (This part has been midified from SPM.py to fit the rest of the script)
    # Find the index of the vertical stringer at x_force
    idx_x_force = np.where(x_Strs_vert == x_force)[0][0]  # Get index

    # Find the top-most y-dof for this vertical stringer
    top_y_dof = y_dof_matrix[np.sum(NrStr_vert_segments_on_row[:idx_x_force+1])-1, 2]  # Last row, middle DOF

    ForceDof = top_y_dof
    f = np.zeros(NrDofs)
    f[ForceDof] = Force_dim

# Process imposed displacements
sR = np.zeros((NrFixedDofs, NrDofs))
fR = np.zeros(NrFixedDofs)
for i in range(NrFixedDofs):
    k = FixedDof[i]
    sR[i, :] = s[k, :]
    fR[i] = f[k]
    s[k, :] = 0.0
    s[k, k] = 1.0
    f[k] = Disp[i]

# Solve system
u = np.linalg.solve(s, f)

# Compute stringer forces
N1 = (StrEA / Strl) * (-4 * u[StrDof[:, 0]] + 6 * u[StrDof[:, 1]] - 2 * u[StrDof[:, 2]])
N2 = (StrEA / Strl) * (2 * u[StrDof[:, 0]] - 6 * u[StrDof[:, 1]] + 4 * u[StrDof[:, 2]])

# Compute panel forces
Taut = PanelGt * ((u[PanelDof[:, 1]] - u[PanelDof[:, 0]]) / Panelb + (u[PanelDof[:, 3]] - u[PanelDof[:, 2]]) / Panela)
Fxy1 = -Panela * Taut
Fxy2 = Panela * Taut
Fyx1 = -Panelb * Taut
Fyx2 = Panelb * Taut

# -- End of the part taken from SPM.py by Hoogenboom --

# Dimensioning of reinforcement
fyd = fyk # [Pa] (No safety factor was used for the purpose of the thesis)

# Tension reinforcement in stringers
print("Reinforcement area needed in stringers")
N_max_list = np.zeros(NrStrs)
A_s_list = np.zeros(NrStrs)
for i in range(NrStrs):  
    N_max = max(N1[i], N2[i])  
    if N_max <= 0:
        A_s = 0  
    else:
        A_s = (N_max/fyd)*(1000000)  

    N_max_list[i] = N_max
    if A_s < (2*(8/2)**2*np.pi):
        A_s = 0
    A_s_list[i] = A_s # [mm^2]
print(A_s_list) # Reinforcement area needed for each tension stringer

# Distributed reinforcement in panels
print("Reinforcement area needed in panels")
A_s_distr_list = np.zeros(NrPanels)
for i in range(NrPanels):
    rho = abs(Taut[i])/fyk # Steel ratio 
    A_s_distr = rho*t*1000000 # [mm2/m]
    if A_s_distr < 400:
        A_s_distr = 400 # [mm2/m]
    A_s_distr_list[i] = A_s_distr 
print(A_s_distr_list) # Distributed reinforcement needed for each panel

# Calculate reinforcement volume
V_stringer_reinf = A_s_list * Strl*1000  # [mm³]
V_panel_reinf = A_s_distr_list*1e-6 * Panela * Panelb # [m3] 
if hole == "yes": 
    for i in range(NrPanels):
        if (i == 1):
            V_panel_reinf[i] += c_edge*Panela[i]*2*400*1e-6
        elif (i == 0) or (i == 2) or (i == 5) or (i == 8):
            V_panel_reinf[i] += (c_supports*Panelb[i] + c_edge*Panela[i])*400*1e-6
        elif (i == 3) or (i == 4): 
            V_panel_reinf[i] += (c_edge*Panelb[i] + c_supports*Panelb[i])*400*1e-6
        elif (i == 6) or (i == 7):
            V_panel_reinf[i] += (c_edge*Panela[i]*2)*400*1e-6
else:
    for i in range(NrPanels):
        V_panel_reinf[i] += (c_supports*Panelb[i] + 2*c_edge*Panela[i])*400*1e-6
V_stringer_reinf = V_stringer_reinf*1e-9 # [m3]
V_panel_reinf = V_panel_reinf*4 # [m3]

Stringer_reinf_volume = np.sum(V_stringer_reinf)*1e6  # [10^-6 m3] Total stringer reinforcement volume
Panel_reinf_volume = np.sum(V_panel_reinf)*1e6 # [10^-6 m3] Total distributed reinforcement volume 

Total_reinf_volume = (Stringer_reinf_volume + Panel_reinf_volume) # [10^-6 m3] Total reinforcement volume
print("Total reinforcement volume")
print(Total_reinf_volume)


# Beräkningar för aktuell last
# -- The following calculations are taken from the script SPM.py by Hoogenboom (with some additions) --

# Initialize stiffness matrix
s = np.zeros((NrDofs, NrDofs))

# Assemble stringers
sl = np.array([[4, -6, 2], [-6, 12, -6], [2, -6, 4]])
for i in range(NrStrs):
    c = StrEA[i] / Strl[i]
    for k in range(3):
        kk = StrDof[i, k]
        for l in range(3):
            ll = StrDof[i, l]
            s[kk, ll] += c * sl[k, l]

# Assemble panels
for i in range(NrPanels):
    a = Panela[i] / Panelb[i]
    e = 1 / a
    c = PanelGt[i]
    sl = np.array([[a, -a, 1, -1], [-a, a, -1, 1], [1, -1, e, -e], [-1, 1, -e, e]])
    for k in range(4):
        kk = PanelDof[i, k] 
        for l in range(4):
            ll = PanelDof[i, l] 
            s[kk, ll] += c * sl[k, l]

# Process imposed forces (This part has been midified from SPM.py to fit the rest of the script)
    # Find the index of the vertical stringer at x_force
    idx_x_force = np.where(x_Strs_vert == x_force)[0][0]  # Get index

    # Find the top-most y-dof for this vertical stringer
    top_y_dof = y_dof_matrix[np.sum(NrStr_vert_segments_on_row[:idx_x_force+1])-1, 2]  # Last row, middle DOF

    ForceDof = top_y_dof
    f = np.zeros(NrDofs)
    f[ForceDof] = Force

# Process imposed displacements
sR = np.zeros((NrFixedDofs, NrDofs))
fR = np.zeros(NrFixedDofs)
for i in range(NrFixedDofs):
    k = FixedDof[i]
    sR[i, :] = s[k, :]
    fR[i] = f[k]
    s[k, :] = 0.0
    s[k, k] = 1.0
    f[k] = Disp[i]

# Solve system
u = np.linalg.solve(s, f)
print("dof displacements")
for i in range(len(u)):
    print(f"{i} {u[i]:.5f}")

# Compute stringer forces
print("stringer forces")
N1 = (StrEA / Strl) * (-4 * u[StrDof[:, 0]] + 6 * u[StrDof[:, 1]] - 2 * u[StrDof[:, 2]])
N2 = (StrEA / Strl) * (2 * u[StrDof[:, 0]] - 6 * u[StrDof[:, 1]] + 4 * u[StrDof[:, 2]])
for i in range(NrStrs):
    print(f"{i + 1} {N1[i]:.1f} {N2[i]:.1f}")

# Compute panel forces
print("panels")
Taut = PanelGt * ((u[PanelDof[:, 1]] - u[PanelDof[:, 0]]) / Panelb + (u[PanelDof[:, 3]] - u[PanelDof[:, 2]]) / Panela)
Fxy1 = -Panela * Taut
Fxy2 = Panela * Taut
Fyx1 = -Panelb * Taut
Fyx2 = Panelb * Taut

for i in range(NrPanels):
    print(f"{i + 1} {Taut[i]:10.3f} {Fxy1[i]:10.3f} {Fxy2[i]:7.3f} {Fyx1[i]:7.3f} {Fyx2[i]:7.3f}")


# Compute support reactions
print("support reactions")
Reactions = np.dot(sR, u) - fR  # Vectorized computation

for i in range(NrFixedDofs):
    print(f"{i + 1} {FixedDof[i]:3d} {Reactions[i]:10.4f}")

# -- End of the part taken from SPM.py by Hoogenboom --


# Calculation of stresses in bottom tensile reinforcement
N_max_list = np.zeros(NrStrs)
for i in range(NrStrs):  
    N_max = max(N1[i], N2[i])  
    N_max_list[i] = N_max

if hole == "yes":
    N = N_max_list[1]
    A = A_s_list[1]
    sigma_s = N/A # MPa
else:
    N = N_max_list[0]
    A = A_s_list[0]
    sigma_s = N/A # MPa


# Calculation of stresses in reinforcement
N_1_list = np.zeros(NrStrs)
N_2_list = np.zeros(NrStrs)
sigma_1_list = np.zeros(NrStrs)
sigma_2_list = np.zeros(NrStrs)
for i in range(NrStrs):  
    N_1 = N1[i]
    N_1_list[i] = N_1
    N_2 = N2[i]
    N_2_list[i] = N_2
    N_1 = N_1_list[i]
    A = A_s_list[i]
    if A == 0:
        continue
    sigma_s_1 = N_1/A # [MPa]
    sigma_1_list[i] = sigma_s_1

    N_2 = N_2_list[i]
    A = A_s_list[i]
    sigma_s_2 = N_2/A # [MPa]
    sigma_2_list[i] = sigma_s_2


# PLOT: showing the nodal coordinates 
import matplotlib.pyplot as plt
# Beam outline
fig, ax = plt.subplots()
ax.set_xlim(0, width)
ax.set_ylim(0, height)
# Nodes
ax.scatter(x_nodes_strs, y_nodes_strs, color='red', label="Nodes")
# Stringers 
for i in range(NrStrs_horiz_full):
    for j in range(NrStr_horiz_segments_on_row[i]):
        if hole == "yes" and j >= 2 and i <= 1:
            x1, x2 = x_Strs_vert[j+1], x_Strs_vert[j + 2] 
        elif hole == "yes" and j == 1 and i <= 1:
            x1, x2 = x_Strs_vert[j], x_Strs_vert[j + 2]
        elif hole == "yes" and j > 1 and i <= 1:
            x1, x2 = x_Strs_vert[j+1], x_Strs_vert[j + 2]
        else: 
            x1, x2 = x_Strs_vert[j], x_Strs_vert[j + 1]
        y = y_Strs_horiz[i]
        ax.plot([x1, x2], [y, y], color='blue')
for i in range(NrStrs_vert_full):
    for j in range(NrStr_vert_segments_on_row[i]):
        if hole == "yes" and i == 2:
            y1, y2 = y_Strs_horiz[j+2], y_Strs_horiz[j + 3]
        else:
            y1, y2 = y_Strs_horiz[j], y_Strs_horiz[j + 1]

        x = x_Strs_vert[i]
        ax.plot([x, x], [y1, y2], color='blue')
# Opening
if hole == "yes":
    hole_rect = plt.Rectangle((x_bottomleft_hole, y_bottomleft_hole), width_hole, height_hole,
                              edgecolor='black', facecolor='none', linestyle='dashed', label="Hole")
    ax.add_patch(hole_rect)
# Node coordinates
for i in range(NrNodes):
    ax.text(x_nodes_strs[i], y_nodes_strs[i], f'({x_nodes_strs[i]:.2f}, {y_nodes_strs[i]:.2f})', 
            fontsize=8, ha='right', color='red')
# Force arrow
force_arrow_x = x_force
force_arrow_y = height  # Nod vid toppen av balken
arrow_length = height * 0.15  # Förkortad pil
ax.annotate("", 
            xy=(force_arrow_x, force_arrow_y - 0.001),  # Pilen slutar vid noden
            xytext=(force_arrow_x, force_arrow_y + arrow_length),  # Starta lite ovanför
            arrowprops=dict(facecolor='red', edgecolor='red', linewidth=2, headwidth=8, headlength=8))
ax.text(force_arrow_x, force_arrow_y + arrow_length - 0.1, f'{abs(Force):.0f} N', 
        color='red', fontsize=10, ha='center',
        bbox=dict(facecolor='white', alpha=0.6, edgecolor='none'))

ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
ax.set_title("Stringer Layout with Node Coordinates and Applied Force")
ax.set_xlim(-0.5, width + 0.5)
ax.set_ylim(-0.5, height + 0.5)
ax.set_aspect('equal')
plt.grid(True, linestyle="--", alpha=0.5)
ax.set_xlabel("Width (m)")
ax.set_ylabel("Height (m)")
ax.legend()
ax.set_title("Stringer Panel Method - Geometry")
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
plt.show()


# PLOT: Numbering of stringers and panels
import matplotlib.pyplot as plt
# Beam outline
fig, ax = plt.subplots()
ax.set_xlim(0, width)
ax.set_ylim(0, height)
# Stringers
Str_coord_list = np.zeros((NrStrs, 4))
S_counter = 1  
for i in range(NrStrs_horiz_full):
    for j in range(NrStr_horiz_segments_on_row[i]):
        if hole == "yes" and j >= 2 and i <= 1:
            x1, x2 = x_Strs_vert[j+1], x_Strs_vert[j + 2]
        elif hole == "yes" and j == 1 and i <= 1:
            x1, x2 = x_Strs_vert[j], x_Strs_vert[j + 2]
        elif hole == "yes" and j > 1 and i <= 1:
            x1, x2 = x_Strs_vert[j+1], x_Strs_vert[j + 2]
        else: 
            x1, x2 = x_Strs_vert[j], x_Strs_vert[j + 1]
        y = y_Strs_horiz[i]
        ax.plot([x1, x2], [y, y], color='blue')
        # Add number to stringer segment
        ax.text((x1 + x2) / 2, y, f'S{S_counter}', fontsize=10, ha='center', color='blue')
        Str_coord_list[S_counter-1,0] = x1
        Str_coord_list[S_counter-1,1] = x2
        Str_coord_list[S_counter-1,2] = y
        Str_coord_list[S_counter-1,3] = y
        S_counter += 1
for i in range(NrStrs_vert_full):
    for j in range(NrStr_vert_segments_on_row[i]):
        if hole == "yes" and i == 2:
            y1, y2 = y_Strs_horiz[j+2], y_Strs_horiz[j + 3]
        else:
            y1, y2 = y_Strs_horiz[j], y_Strs_horiz[j + 1]
        x = x_Strs_vert[i]
        ax.plot([x, x], [y1, y2], color='blue')
        # Add number to stringer segment
        ax.text(x, (y1 + y2) / 2, f'S{S_counter}', fontsize=10, ha='center', color='blue')
        Str_coord_list[S_counter-1,0] = x
        Str_coord_list[S_counter-1,1] = x
        Str_coord_list[S_counter-1,2] = y1
        Str_coord_list[S_counter-1,3] = y2
        S_counter += 1
# Opening
if hole == "yes":
    hole_rect = plt.Rectangle((x_bottomleft_hole, y_bottomleft_hole), width_hole, height_hole, edgecolor='black', facecolor='none', linestyle='dashed', label="Hole")
    ax.add_patch(hole_rect)
# Force arrow
force_arrow_x = x_force
force_arrow_y = height  
arrow_length = height * 0.15  
ax.annotate("", xy=(force_arrow_x, force_arrow_y - 0.001), xytext=(force_arrow_x, force_arrow_y + arrow_length), arrowprops=dict(facecolor='red', edgecolor='red', linewidth=2, headwidth=8, headlength=8))
# Add numbers to panels
panel_centers = np.zeros((NrPanels,2))
P_counter = 1  
for i in range(NrStr_vert_segments):  
    for j in range(NrStr_horiz_segments_on_row[i]):  
        if hole == "yes" and j >= 2 and i <= 1:
            x1, x2 = x_Strs_vert[j+1], x_Strs_vert[j + 2] 
        elif hole == "yes" and i == 1 and j == 1: 
            continue
        elif hole == "yes" and j == 1 and i <= 1:
            x1, x2 = x_Strs_vert[j], x_Strs_vert[j + 2]
        elif hole == "yes" and j > 1 and i <= 1:
            x1, x2 = x_Strs_vert[j+1], x_Strs_vert[j + 2]
        else: 
            x1, x2 = x_Strs_vert[j], x_Strs_vert[j + 1]
        y1 = y_Strs_horiz[i]
        y2 = y_Strs_horiz[i + 1]
        # Add number to the panel
        panel_centers[(P_counter-1)] = (x1 + x2) / 2, (y1 + y2) / 2
        ax.text((x1 + x2) / 2, (y1 + y2) / 2, f'P{P_counter}', fontsize=10, ha='center', color='green')
        P_counter += 1

ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
ax.set_title("Stringer Layout with Node Coordinates and Applied Force")
ax.set_xlim(-0.5, width + 0.5)
ax.set_ylim(-0.5, height + 0.5)
ax.set_aspect('equal')
plt.grid(True, linestyle="--", alpha=0.5)
ax.set_xlabel("Width (m)")
ax.set_ylabel("Height (m)")
ax.set_title("Stringer Panel Method - Geometry")
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
plt.show()

# PLOT: Normal forces in stringers and shear forces in panels
import numpy as np
import matplotlib.pyplot as plt
def plot_lines_with_forces(Str_coord_list, N1, N2, Taut, panel_centers, scale=3e-7, shear_scale=0.1, arrow_offset=0.08):
    fig, ax = plt.subplots()
    for i in range(len(Str_coord_list)):
        x1, x2, y1, y2 = Str_coord_list[i]
        dx = x2 - x1
        dy = y2 - y1
        length = np.sqrt(dx**2 + dy**2)
        if length == 0:
            continue  
        normal_x = -dy / length  
        normal_y = dx / length
        # Scale the normal forces for visualization
        fx1 = normal_x * N1[i] * scale
        fy1 = normal_y * N1[i] * scale
        fx2 = normal_x * N2[i] * scale
        fy2 = normal_y * N2[i] * scale
        ax.plot([x1, x2], [y1, y2], 'b-', linewidth=2) # Structural line
        fill_x = [x1, x1 + fx1, x2 + fx2, x2]
        fill_y = [y1, y1 + fy1, y2 + fy2, y2]
        ax.fill(fill_x, fill_y, color='lightgray', alpha=0.5)
    # Shear stress in panels
    for i, (px, py) in enumerate(panel_centers):
        tau = Taut[i]
        shear_dir = np.sign(tau)
        s = shear_scale
        square_x = [px - s, px + s, px + s, px - s, px - s]
        square_y = [py - s, py - s, py + s, py + s, py - s]
        ax.plot(square_x, square_y, 'k-', linewidth=1)
        # Define arrows around the square 
        arrows = [
            (px, py + s + arrow_offset, -shear_dir * s, 0),  # Top (left)
            (px - s - arrow_offset, py, 0, shear_dir * s),   # Left (up)
            (px, py - s - arrow_offset, shear_dir * s, 0),  # Bottom (right)
            (px + s + arrow_offset, py, 0, -shear_dir * s)  # Right (down)
        ]
        for x, y, dx, dy in arrows:
            ax.arrow(x, y, dx, dy, head_width=s/2, head_length=s/2, fc='k', ec='k')
        # Display shear stress value in the center of the square
        ax.text(px, py, f'{(tau*1e-3):.0f}', fontsize=8, ha='center', va='center', color='black', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_title("Normal Force Distribution")
    ax.set_aspect('equal')
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.show()
plot_lines_with_forces(Str_coord_list, N1, N2, Taut, panel_centers)


# PLOT: Stresses in reinforcement
import numpy as np
import matplotlib.pyplot as plt
def plot_lines_with_stresses(Str_coord_list, sigma_1_list, sigma_2_list, A_s_list, width, height,x_bottomleft_hole,y_bottomleft_hole,width_hole,height_hole, scale=5e-4):
    fig, ax = plt.subplots()
    S_num = 0
    for i in range(len(Str_coord_list)):
        S_num +=1
        if A_s_list[i] == 0:
            continue  
        x1, x2, y1, y2 = Str_coord_list[i]
        dx = x2 - x1
        dy = y2 - y1
        length = np.sqrt(dx**2 + dy**2)
        if length == 0:
            continue  
        normal_x = -dy / length  
        normal_y = dx / length
        # Scale the stress values for visualization
        fx1 = normal_x * sigma_1_list[i] * scale
        fy1 = normal_y * sigma_1_list[i] * scale
        fx2 = normal_x * sigma_2_list[i] * scale
        fy2 = normal_y * sigma_2_list[i] * scale
        ax.plot([x1, x2], [y1, y2], 'k-', linewidth=2) # Structural line
        fill_x = [x1, x1 + fx1, x2 + fx2, x2]
        fill_y = [y1, y1 + fy1, y2 + fy2, y2]
        ax.fill(fill_x, fill_y, color='lightblue', alpha=0.5)
        mid_x = (x1 + x2) / 2
        mid_y = (y1 + y2) / 2
        ax.text(mid_x, mid_y - 0.15, f'$S_{S_num}$', fontsize=10, ha='center', color='black') # Add stringer number 
        # Beam outline 
        beam_x = [0, width, width, 0, 0]
        beam_y = [0, 0, height, height, 0]
        ax.plot(beam_x, beam_y, color='gray', linestyle='--', linewidth=1)
    ax.set_xlabel("[m]")
    ax.set_ylabel("[m]")
    ax.set_title("Stress Distribution in Reinforcement")
    ax.set_aspect('equal')
    plt.show()
plot_lines_with_stresses(Str_coord_list, sigma_1_list, sigma_2_list, A_s_list, width, height,x_bottomleft_hole,y_bottomleft_hole,width_hole,height_hole)
