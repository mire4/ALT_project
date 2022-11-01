import numpy as np

def levenshtein_matriz(x, y, threshold=None):
    # esta versión no utiliza threshold, se pone porque se puede
    # invocar con él, en cuyo caso se ignora
    lenX, lenY = len(x), len(y)
    D = np.zeros((lenX + 1, lenY + 1), dtype=np.int)
    for i in range(1, lenX + 1): # Rellena la primera fila
        D[i][0] = D[i - 1][0] + 1
    for j in range(1, lenY + 1): # Para cada columna
        D[0][j] = D[0][j - 1] + 1 # Rellenas su vertical
        for i in range(1, lenX + 1): # Para cada fila
            D[i][j] = min(
                D[i - 1][j] + 1,
                D[i][j - 1] + 1,
                D[i - 1][j - 1] + (x[i - 1] != y[j - 1]),
            )
    return D[lenX, lenY]

def levenshtein_edicion(x, y, threshold=None):
    # a partir de la versión levenshtein_matriz
    lenX, lenY = len(x), len(y)
    D = np.zeros((lenX + 1, lenY + 1), dtype=np.int)
    for i in range(1, lenX + 1): # Rellena la primera fila
        D[i][0] = D[i - 1][0] + 1
    for j in range(1, lenY + 1): # Para cada columna
        D[0][j] = D[0][j - 1] + 1 # Rellenas su vertical
        for i in range(1, lenX + 1): # Para cada fila
            D[i][j] = min(
                D[i - 1][j] + 1,
                D[i][j - 1] + 1,
                D[i - 1][j - 1] + (x[i - 1] != y[j - 1]),
            )
    indX, indY = D.shape[0] - 1, D.shape[1] - 1
    op = []
    while indX > 0 and indY > 0:
        ins = D[indX, indY - 1]
        bor = D[indX - 1, indY]
        sus = D[indX - 1, indY - 1]
        opMin = min(ins, bor, sus)
        if (sus == opMin):
            op.append('sus')
            indX -= 1
            indY -= 1
        elif (ins == opMin):
            op.append('ins')
            indY -= 1
        elif (bor == opMin):
            op.append('bor')
            indX -= 1
    # Caso llego a una pared de la matriz
    while indY > 0:
        op.append('ins')
        indY -= 1
    while indX > 0:
        op.append('bor')
        indX -= 1

    # Recuperación del camino
    op = op[::-1]
    indX, indY = 0, 0
    secOp = []
    for o in op:
        if (o == 'ins'):
            secOp.append(('', y[indY]))
            indY += 1
        elif (o == 'bor'):
            secOp.append((x[indX], ''))
            indX += 1
        elif (o == 'sus'):
            secOp.append((x[indX], y[indY]))
            indX += 1
            indY += 1
    return D[lenX, lenY], secOp

def levenshtein_reduccion(x, y, threshold=None):
    lenX, lenY = len(x), len(y)
    vcurrent = np.zeros(lenX, dtype=int)
    vprev = np.arange(1, lenX + 1, dtype=int)
    for j in range(lenY):
        for i in range(lenX):
            if (i == 0):
                vcurrent[0] = min(
                    j + (x[i] != y[j]),
                    vprev[i] + 1
                )
            else:
                vcurrent[i] = min(
                    vcurrent[i - 1] + 1,
                    vprev[i] + 1,
                    vprev[i - 1] + (x[i] != y[j])
                )
        vprev, vcurrent = vcurrent, vprev
    return vprev[-1]
                
def levenshtein(x, y, threshold):
    lenX, lenY = len(x), len(y)
    difX = lenX - lenY
    vcurrent = np.zeros(lenX, dtype=int)
    vprev = np.arange(1, lenX + 1, dtype=int)
    for j in range(lenY):
        if (difX + j >= 0 and vprev[difX + j] > threshold):
            return threshold + 1
        for i in range(lenX):
            if (i == 0):
                vcurrent[0] = min(
                    j + (x[i] != y[j]),
                    vprev[i] + 1
                )
            else:
                vcurrent[i] = min(
                    vcurrent[i - 1] + 1,
                    vprev[i] + 1,
                    vprev[i - 1] + (x[i] != y[j])
                )
        vprev, vcurrent = vcurrent, vprev
    return vprev[-1]

def levenshtein_cota_optimista(x, y, threshold):
    return 0 # COMPLETAR Y REEMPLAZAR ESTA PARTE

def damerau_restricted_matriz(x, y, threshold=None):
    # completar versión Damerau-Levenstein restringida con matriz
    lenX, lenY = len(x), len(y)
    # COMPLETAR
    return D[lenX, lenY]

def damerau_restricted_edicion(x, y, threshold=None):
    # partiendo de damerau_restricted_matriz añadir recuperar
    # secuencia de operaciones de edición
    return 0,[] # COMPLETAR Y REEMPLAZAR ESTA PARTE

def damerau_restricted(x, y, threshold=None):
    # versión con reducción coste espacial y parada por threshold
     return min(0,threshold+1) # COMPLETAR Y REEMPLAZAR ESTA PARTE

def damerau_intermediate_matriz(x, y, threshold=None):
    # completar versión Damerau-Levenstein intermedia con matriz
    return D[lenX, lenY]

def damerau_intermediate_edicion(x, y, threshold=None):
    # partiendo de matrix_intermediate_damerau añadir recuperar
    # secuencia de operaciones de edición
    # completar versión Damerau-Levenstein intermedia con matriz
    return 0,[] # COMPLETAR Y REEMPLAZAR ESTA PARTE
    
def damerau_intermediate(x, y, threshold=None):
    # versión con reducción coste espacial y parada por threshold
    return min(0,threshold+1) # COMPLETAR Y REEMPLAZAR ESTA PARTE

opcionesSpell = {
    'levenshtein_m': levenshtein_matriz,
    'levenshtein_r': levenshtein_reduccion,
    'levenshtein':   levenshtein,
    'levenshtein_o': levenshtein_cota_optimista,
    'damerau_rm':    damerau_restricted_matriz,
    'damerau_r':     damerau_restricted,
    'damerau_im':    damerau_intermediate_matriz,
    'damerau_i':     damerau_intermediate
}

opcionesEdicion = {
    'levenshtein': levenshtein_edicion,
    'damerau_r':   damerau_restricted_edicion,
    'damerau_i':   damerau_intermediate_edicion
}

