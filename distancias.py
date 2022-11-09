import numpy as np

def levenshtein_matriz(x, y, threshold=None):
    """
    Cálculo distancia de levenshtein con matriz
    """
    lenX, lenY = len(x), len(y)
    D = np.zeros((lenX + 1, lenY + 1), dtype=np.int)
    for i in range(1, lenX + 1): # Rellena la primera fila
        D[i][0] = D[i - 1][0] + 1
    for j in range(1, lenY + 1): # Para todas las filas
        D[0][j] = D[0][j - 1] + 1 # Rellenas el primer elemento
        for i in range(1, lenX + 1): # Para el resto de columnas
            D[i][j] = min(
                D[i - 1][j] + 1, # Inserción
                D[i][j - 1] + 1, # Borrado
                D[i - 1][j - 1] + (x[i - 1] != y[j - 1]), # Sustitución
            )
    return D[lenX, lenY]

def levenshtein_edicion(x, y, threshold=None):
    """
    Cálculo distancia de levenshtein con matriz y 
    recorrido
    """
    lenX, lenY = len(x), len(y)
    D = np.zeros((lenX + 1, lenY + 1), dtype=np.int)
    for i in range(1, lenX + 1): # Rellena la primera fila
        D[i][0] = D[i - 1][0] + 1
    for j in range(1, lenY + 1): # Para todas las filas
        D[0][j] = D[0][j - 1] + 1 # Rellenas el primer elemento
        for i in range(1, lenX + 1): # Para el resto de columnas
            D[i][j] = min(
                D[i - 1][j] + 1, # Inserción
                D[i][j - 1] + 1, # Borrado
                D[i - 1][j - 1] + (x[i - 1] != y[j - 1]), # Sustitución
            )

    indX, indY = D.shape[0] - 1, D.shape[1] - 1  # Desde el último elemento
    secOp = []
    while indX > 0 and indY > 0:
        ins = D[indX, indY - 1]
        bor = D[indX - 1, indY]
        sus = D[indX - 1, indY - 1]
        op, decX, decY = elegir_operacion(x, y, indX, indY, ins, bor, sus)
        secOp.append(op)
        indX -= decX
        indY -= decY
    while indY > 0:
        secOp.append(('', y[indY - 1]))
        indY -= 1
    while indX > 0:
        secOp.append((x[indX - 1], ''))
        indX -= 1
    secOp = secOp[::-1]
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
        if (difX + j > 0 and vprev[difX + j - 1] > threshold):
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
    contDict = dict()
    distancePos, distanceNeg = 0, 0
    for letter in x:
        contDict[letter] = contDict.get(letter, 0) + 1
    for letter in y:
        contDict[letter] = contDict.get(letter, 0) - 1
    for value in contDict.values():
        if value > 0:
            distancePos += value
        else:
            distanceNeg += abs(value)
        if distancePos > threshold or distanceNeg > threshold:
            return threshold + 1
    return levenshtein(x, y, threshold)

def damerau_restricted_matriz(x, y, threshold=None):
    # completar versión Damerau-Levenstein restringida con matriz
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
            if (i > 1 and j > 1)\
                and (x[i - 2] == y[j - 1] and x[i - 1] == y[j - 2]):
                D[i][j] = min(
                    D[i][j],
                    D[i - 2][j - 2] + 1
                )
                
    return D[lenX, lenY]

def damerau_restricted_edicion(x, y, threshold=None):
    # partiendo de damerau_restricted_matriz añadir recuperar
    # secuencia de operaciones de edición
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
            if (i > 1 and j > 1)\
                and (x[i - 2] == y[j - 1] and x[i - 1] == y[j - 2]):
                D[i][j] = min(
                    D[i][j],
                    D[i - 2][j - 2] + 1
                )
    indX, indY = D.shape[0] - 1, D.shape[1] - 1
    secOp = []
    while indX > 0 and indY > 0:
        ins = D[indX, indY - 1]
        bor = D[indX - 1, indY]
        sus = D[indX - 1, indY - 1]
        if (indX > 1 and indY > 1)\
                and (x[indX - 2] == y[indY - 1] and x[indX - 1] == y[indY - 2]):
            int = D[indX - 2, indY - 2]
        else:
            int = float('inf')
        op, decX, decY = elegir_operacion(x, y, indX, indY, ins, bor, sus, int)
        secOp.append(op)
        indX -= decX
        indY -= decY
    while indY > 0:
        secOp.append(('', y[indY - 1]))
        indY -= 1
    while indX > 0:
        secOp.append((x[indX - 1], ''))
        indX -= 1
    secOp = secOp[::-1]
    return D[lenX, lenY], secOp

def damerau_restricted(x, y, threshold=None):
    # versión con reducción coste espacial y parada por threshold
    lenX, lenY = len(x), len(y)
    difX = lenX - lenY
    vcurrent = np.zeros(lenX, dtype=int)
    vprev = np.arange(1, lenX + 1, dtype=int)
    vprevprev = 2147483647 * np.ones(lenX, dtype=int) # inf
    for j in range(lenY):
        for i in range(lenX):
            if (difX + j > 0 and vprev[difX + j - 1] > threshold):
                return threshold + 1
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
            if ((i > 0 and j > 0)\
                and (x[i - 1] == y[j] and x[i] == y[j - 1])):
                if (i == 1):
                    vcurrent[i] = min(
                        vcurrent[i],
                        j               # j - 1 + 1
                    )
                else:
                    vcurrent[i] = min(
                        vcurrent[i],
                        vprevprev[i - 2] + 1
                    )
        vprevprev, vprev = vprev, vprevprev
        vprev, vcurrent = vcurrent, vprev
    return vprev[-1]

def damerau_intermediate_matriz(x, y, threshold=None):
    # completar versión Damerau-Levenstein intermedia con matriz
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
            if (i > 1 and j > 1)\
                and (x[i - 2] == y[j - 1] and x[i - 1] == y[j - 2]):
                D[i][j] = min(
                    D[i][j],
                    D[i - 2][j - 2] + 1
                )

            # acb -> ba
            if (i > 2 and j > 1)\
                and (x[i - 3] == y[j - 1] and x[i - 1] == y[j - 2]):
                D[i][j] = min(
                    D[i][j],
                    D[i - 3][j - 2] + 2
                )

            # ba -> acb
            if (i > 1 and j > 2)\
                and (x[i - 2] == y[j - 1] and x[i - 1] == y[j - 3]):
                D[i][j] = min(
                    D[i][j],
                    D[i - 2][j - 3] + 2
                )
                
    return D[lenX, lenY]

def damerau_intermediate_edicion(x, y, threshold=None):
    # partiendo de matrix_intermediate_damerau añadir recuperar
    # secuencia de operaciones de edición
    # completar versión Damerau-Levenstein intermedia con matriz
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
            if (i > 1 and j > 1)\
                and (x[i - 2] == y[j - 1] and x[i - 1] == y[j - 2]):
                D[i][j] = min(
                    D[i][j],
                    D[i - 2][j - 2] + 1
                )
            #acb -> ba
            if (i > 2 and j > 1)\
                and (x[i - 3] == y[j - 1] and x[i - 1] == y[j - 2]):
                D[i][j] = min(
                    D[i][j],
                    D[i - 3][j - 2] + 2
                )

            #ba -> acb
            if (i > 1 and j > 2)\
                and (x[i - 2] == y[j - 1] and x[i - 1] == y[j - 3]):
                D[i][j] = min(
                    D[i][j],
                    D[i - 2][j - 3] + 2
                )

    indX, indY = D.shape[0] - 1, D.shape[1] - 1
    secOp = []
    while indX > 0 and indY > 0:
        ins = D[indX, indY - 1]
        bor = D[indX - 1, indY]
        sus = D[indX - 1, indY - 1]
        if (indX > 1 and indY > 1)\
                and (x[indX - 2] == y[indY - 1] and x[indX - 1] == y[indY - 2]):
            int = D[indX - 2, indY - 2]
        else:
            int = float('inf')
        if (indX > 2 and indY > 1)\
                and (x[indX - 3] == y[indY - 1] and x[indX - 1] == y[indY - 2]):
            int3_2 = D[indX - 3][indY - 2]
        else:
            int3_2 = float('inf')
        if (indX > 1 and indY > 2)\
                and (x[indX - 2] == y[indY - 1] and x[indX - 1] == y[indY - 3]):
            int2_3 = D[indX - 2][indY - 3]
        else:
            int2_3 = float('inf')
        op, decX, decY = elegir_operacion(x, y, indX, indY, ins, bor, sus, int, int3_2, int2_3)
        secOp.append(op)
        indX -= decX
        indY -= decY
    while indY > 0:
        secOp.append(('', y[indY - 1]))
        indY -= 1
    while indX > 0:
        secOp.append((x[indX - 1], ''))
        indX -= 1
    secOp = secOp[::-1]
    return D[lenX, lenY], secOp
    
def damerau_intermediate(x, y, threshold=None):
    # versión con reducción coste espacial y parada por threshold
    lenX, lenY = len(x), len(y)
    difX = lenX - lenY
    vcurrent = np.zeros(lenX, dtype=int)
    vprev = np.arange(1, lenX + 1, dtype=int)
    vprevprev = 2147483647 * np.ones(lenX, dtype=int) # inf
    vprevprevprev = 2147483647 * np.ones(lenX, dtype=int) # inf
    for j in range(lenY):
        for i in range(lenX):
            if (difX + j > 0 and vprev[difX + j - 1] > threshold):
                return threshold + 1
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
            if ((i > 0 and j > 0)\
                and (x[i - 1] == y[j] and x[i] == y[j - 1])):
                if (i == 1):
                    vcurrent[i] = min(
                        vcurrent[i],
                        j               # j - 1 + 1
                    )
                else:
                    vcurrent[i] = min(
                        vcurrent[i],
                        vprevprev[i - 2] + 1
                    )
            #acb -> ba
            if (i > 1 and j > 0)\
                and (x[i - 2] == y[j] and x[i] == y[j - 1]):
                if (i == 2):
                    vcurrent[i] = min(
                        vcurrent[i],
                        j + 1           # j - 1 + 2
                    )
                else:
                    vcurrent[i] = min(
                        vcurrent[i],
                        vprevprev[i - 3] + 2
                    )
            #ba -> acb
            if (i > 0 and j > 1)\
                and (x[i - 1] == y[j] and x[i] == y[j - 2]):
                if (i == 1):
                    vcurrent[i] = min(
                        vcurrent[i],
                        j               # j - 2 + 2
                    )
                else:
                    vcurrent[i] = min(
                        vcurrent[i],
                        vprevprevprev[i - 2] + 2
                    )
        vprevprevprev, vprevprev = vprevprev, vprevprevprev
        vprevprev, vprev = vprev, vprevprev
        vprev, vcurrent = vcurrent, vprev
    return vprev[-1]

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

def elegir_operacion(x, y, indX, indY, dIns, dBor, dSus, dInt=float('inf'), dInt3_2=float('inf'), dInt2_3=float('inf')):
    opMin = min(dIns, dBor, dSus, dInt, dInt3_2, dInt2_3)
    if (dInt3_2 == opMin):
        return (str(x[indX - 3:indX]), str(y[indY - 2:indY])), 3, 2
    elif (dInt2_3 == opMin):
        return (str(x[indX - 2:indX]), str(y[indY - 3:indY])), 2, 3
    elif (dInt == opMin):
        return (str(x[indX - 2:indX]), str(y[indY - 2:indY])), 2, 2
    elif (dSus == opMin):
        return (x[indX - 1], y[indY - 1]), 1, 1
    elif (dIns == opMin):
        return ('', y[indY - 1]), 0, 1
    elif (dBor == opMin):
        return (x[indX - 1], ''), 1, 0
    else:
        print("Error en edición")
        exit()