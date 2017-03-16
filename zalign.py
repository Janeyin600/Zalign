#!/usr/bin/env python2

#*********************
#      zalign.py      
# written by Jane Yin 
#*********************

import sys
import os
import math
import string
import re

def align_exception():
  print ''
  print '********************************************************************************************'
  print '       Welcome to Zalign!'
  print 'Usage: Please provide the name of the PDB file (-f) and the indices of three receptor atoms'
  print '       that you think should lie in the same XY plane after the Z-axis alignment (-a1, -a2, -a3).'
  print '       Pick up a fourth receptor atom that is supposed to have more negative Z coordinates'
  print '       after alignment (-a4), so that the molecule will point to the desired +Z or -Z direction.'
  print '       In addition, specify the rotation angle (in degrees) around the z axis (-r) or leave it as 0.'
  print '       The new origin will be read from the APR input file (the G1 entry).'
  print '\n'           
  print 'For example: python zalign.py -f pdb/oa_cba.pdb -a1 :OCT@O3 -a2 :OCT@O7 -a3 :OCT@O8 -a4 :OCT@O11 -r 85 -i ../apr.in'
  print '********************************************************************************************'
  print '\n'
  return

# Check if the user input is an integer, a float number or a string
def ismyinstance(val_type, val, str):
  # if no input
  if not val: 
        return 0
  if val_type == 'float': 
    try:
        float(val)
    except ValueError:
        print  val
        print ('Please enter a float value for %s.'%(str))
        sys.exit()
  elif val_type == 'int':
    try:
        int(val)
    except ValueError:
        return False
        sys.exit()
  elif (val_type == 'string'):
        return val
  
  if (float(val) < 0):
        print  val
        print ('Please enter a non-negative value for %s.'%(str))
        sys.exit()
   
  if (val_type == 'float'):
        return float(val)
  elif (val_type == 'int'):
        return int(val) 

def find_index(usr_input, pdbfile):
    """
    Find the atom index of an atom in a PDB file.
    :param name: atom name (in PDB file)
    :param resname: atom resname (in PDB file)
    :return:
    """
    atom = usr_input.split('@')[1]
    residue = usr_input.split('@')[0][1:]
    resname  = 'None'
    res_id = 99999

    if ismyinstance('int',residue,'N/A'):
      # residue number  was provided
      res_id = int(residue)
    else:
      resname = residue

    flag = 0 

    with open(pdbfile) as pdb_file:
        lines = (line.rstrip('\n') for line in pdb_file)
        lines = list(line for line in pdb_file) 
 
    for i in range(0, len(lines)):
        newline = lines[i].split()
        if newline[0] == 'ATOM' or newline[0] == 'HETATM':
            if (newline[3] == resname or int(lines[i][22:27].strip())==res_id) and newline[2] == atom:
                flag = 1
                break
    pdb_file.close()

    if not flag :
        print ('%s cannot be found.'%(usr_input))
        sys.exit()    

    return int(newline[1])


# Check if a file exists
def check_file(pdb_file):
  try:
    f = open(pdb_file,'r');
  except IOError:
    print 'Error: File does not exist'
    sys.exit()
  f.close()  

#release the memory of a list
def release_list(a):
  del a[:]
  del a

################################
#      zalign starts          #
                          
if (len(sys.argv)!=15):
   align_exception()
   print 'Aborted. Check to see if you missed any flags.'
   sys.exit()

#read arguments from the command line

arg_list = [1,3,5,7,9,11,13]
for i in arg_list:
  if (sys.argv[i] == '-f'):
    pdb_file = sys.argv[i+1]
  elif (sys.argv[i] == '-a1'):
    H1 = sys.argv[i+1]
  elif (sys.argv[i] == '-a2'):
    H2 = sys.argv[i+1]
  elif (sys.argv[i] == '-a3'):
    H3 = sys.argv[i+1]
  elif (sys.argv[i] == '-a4'):
    H4 = sys.argv[i+1]
  elif (sys.argv[i] == '-i'):
    input_file = sys.argv[i+1]
  elif (sys.argv[i] == '-r'):
    rotangle = sys.argv[i+1]  
  else:
    print 'Wrong flags!! Please only use -f, -a1, -a2, -a3, -a4, -r , or -i.'  

ismyinstance('float',rotangle,'-r')

cutoff = 0.05

check_file(pdb_file)
check_file(input_file) 

# look for G1
with open(input_file) as f_in:
   lines = (line.rstrip() for line in f_in)
   lines = list(line for line in lines if line)  # Non-blank lines in a list

flag = 0

for i in range(0, len(lines)):
    lines[i] = re.sub('[\s]', '', lines[i])
    splitline = lines[i].split("=")
    if splitline[0] == 'G1':
	G1 = splitline[1]
	flag = 1

if flag == 0:
    print ('G1 entry cannot be found in %s.\n'%(input_file))  
    sys.exit(1)

idx1 = find_index(H1,pdb_file)
idx2 = find_index(H2,pdb_file)
idx3 = find_index(H3,pdb_file)
idx4 = find_index(H4,pdb_file)
idx_origin = find_index(G1,pdb_file)
  

# read coordinates and other information from the PDB file

resname_list = []
total_atom  = 0
coords = []
atom_namelist = []
resid_list = []
ter_list = []
chain_list = []
head_list = []
col_list = []
col2_list = []
col3_list = []

newPDB_file = open('align_z.pdb', 'w')

with open(pdb_file) as f_in:
  lines = (line.rstrip() for line in f_in)
  lines = list(line for line in lines if line) # Non-blank lines in a list   

for i in range(len(lines)):
  splitdata = lines[i].split()
  # skip the header lines and seperating lines  
  if (splitdata[0]=='ATOM')or(splitdata[0]=='HETATM'):
          total_atom += 1
          coords.append((float(lines[i][30:38].strip()), float(lines[i][38:46].strip()), float(lines[i][46:54].strip())))
          head_list.append(splitdata[0]) 
          atom_namelist.append(lines[i][12:16])             
          resname_list.append(splitdata[3])
          resid_list.append(int(lines[i][22:27].strip()))
          chain_list.append(lines[i][20:22].strip())
          col_list.append(lines[i][55:59].strip())
          col2_list.append(lines[i][60:66].strip())
          col3_list.append(lines[i][74:78].strip()) 
  elif splitdata[0] == 'TER':
    # keep track of the end of a residue
    ter_list.append(total_atom)

flag = 0
print 'Start searching'
# rotate around the X-axis, 1 degree each time 
for i in range (0, 360):
  coords_new = []
  dx = i*math.pi/180.0
  for num1 in range (0, total_atom):
    tmp1 = coords[num1][0]
    tmp2 = math.cos(dx)*coords[num1][1] + (-1)*math.sin(dx)*coords[num1][2]
    tmp3 = math.sin(dx)*coords[num1][1] + math.cos(dx)*coords[num1][2]
    coords_new.append((tmp1,tmp2,tmp3))
     
  #rotate around the Y-axis, 1 degree each time
 
  for j in range (0, 360):
    coords_new2 =[]
    sys.stdout.write('\rScanning... '+ str((i*360+j)/(36*36)) + '%')
    dy = j*math.pi/180.0
    for num2 in range (0, total_atom):
      tmp1 = math.cos(dy)*coords_new[num2][0] + math.sin(dy)*coords_new[num2][2] 
      tmp2 = coords_new[num2][1]
      tmp3 = (-1)*math.sin(dy)*coords_new[num2][0] + math.cos(dy)*coords_new[num2][2]
      coords_new2.append((tmp1,tmp2,tmp3))
    average = (float(coords_new2[idx1-1][2]) + float(coords_new2[idx2-1][2]) + float(coords_new2[idx3-1][2]))/3.0

    diff1 = float(coords_new2[idx1-1][2])-average
    diff2 = float(coords_new2[idx2-1][2])-average
    diff3 = float(coords_new2[idx3-1][2])-average
    if (abs(diff1) < cutoff)and(abs(diff2)<cutoff)and(abs(diff3)<cutoff)and(float(coords_new2[idx4-1][2])<average):
         flag = 1
         print '\n'
         print 'Solution found.'
         break     #break the inner loop 
  
    release_list(coords_new2)

  if (flag==1):
    break   #break the outer loop 
  release_list(coords_new)

if flag == 0:
   print '\nCannot find a solution. Please try a larger cutoff value. (Modify the value of cutoff variable in the script)'
   sys.exit(0) 

# Translate the coordinates according to the new origin
coords_new = []
for i in range(0, total_atom):
  tmpx = float(coords_new2[i][0]) - float(coords_new2[idx_origin-1][0])
  tmpy = float(coords_new2[i][1]) - float(coords_new2[idx_origin-1][1])
  tmpz = float(coords_new2[i][2]) - float(coords_new2[idx_origin-1][2])
  coords_new.append((tmpx, tmpy, tmpz))

# Rotate around the z axis
coords_new3 = []
dz = float(rotangle)*math.pi/180.0

for i in range(0, total_atom):
  tmpx = math.cos(dz)*float(coords_new[i][0]) + (-1)*math.sin(dz)*float(coords_new[i][1]);
  tmpy = math.sin(dz)*float(coords_new[i][0]) + math.cos(dz)*float(coords_new[i][1]);  
  tmpz = float(coords_new[i][2]); #  z coordinate doesn't change
  coords_new3.append((tmpx, tmpy, tmpz))

# Adjust the residue IDs if they do not start with 4
diff =  resid_list[0] - 4
resid_list = [ id - diff for id in resid_list]

# Write the new pdb file
if (flag == 1):
  # Positions for the dummy atoms 
  newPDB_file.write('ATOM      1  Pb  DUM     1       0.000   0.000  -6.000  1.00  0.00\n')
  newPDB_file.write('TER\n')
  newPDB_file.write('ATOM      2  Pb  DUM     2       0.000   0.000 -11.000  1.00  0.00\n')
  newPDB_file.write('TER\n')
  newPDB_file.write('ATOM      3  Pb  DUM     3       0.000   3.500 -14.500  1.00  0.00\n')
  newPDB_file.write('TER\n')
  

  for i in range(total_atom):
      if not chain_list[i]:
         chain_list[i] = ' '
      newPDB_file.write('%-6s%5s %4s %-3s %s%4s'%(head_list[i], i+4, atom_namelist[i],resname_list[i], ' ', resid_list[i]))
      newPDB_file.write('%12.3f%8.3f%8.3f'%(float(coords_new3[i][0]), float(coords_new3[i][1]), float(coords_new3[i][2])))
      newPDB_file.write('%6.2f%6.2f\n'%(float(col_list[i]), float(col2_list[i])))
      #newPDB_file.write('%12s\n'%(col3_list[i]))
      if i+1 in ter_list:      
         newPDB_file.write('%-6s%5s %4s %-3s %s%4s\n'%('TER', i+4, ' ',resname_list[i], ' ', resid_list[i]))

  if total_atom not in ter_list:
     newPDB_file.write('TER\n') 

  newPDB_file.write('END\n')
  print 'A new pdb file, align_z.pdb has been generated.'

f_in.close()
newPDB_file.close() 


