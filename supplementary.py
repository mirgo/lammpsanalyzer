def readlammpstrjpositions(numberofatoms, fileloc): # get matrix with width 3, for x,y,z
  f = open(fileloc, 'r')
  b=0
  positions = []
  for line in tqdm.tqdm(f):
    if b>=9:  # start returning values when passing the header info
      v = ([m.start() for m in re.finditer(' ', line)])[:4] # indices of all spaces in a string
      x = float(line[v[0]:v[1]])
      y = float(line[v[1]:v[2]])
      z = float(line[v[2]:v[3]])
      positions.append([x,y,z])
    b+=1
    if b== 9 + numberofatoms: # reset counter every number of atoms cycle
      b=0
  f.close()
  return np.array(positions)
