import numpy as np
matrix = [
    [394, 378, 397], 
    [388, 376, 427], 
    [418, 425, 397]]
matrix = np.array(matrix)
R,C = np.unravel_index(matrix.argmax(), matrix.shape)
print(R,C, matrix[R][C])
nb_list = [
    [(R-1)%3,C],[(R+1)%3,C],
    [R,(C-1)%3],[R,(C+1)%3]]
nb_val = [matrix[x][y] for x,y in nb_list]
mR,mC = nb_list[nb_val.index(min(nb_val))]
print(mR, mC, matrix[mR][mC])
