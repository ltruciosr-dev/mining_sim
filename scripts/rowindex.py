class caractModel:
    """
    Se define las caracteristicas basica de un modelo de bloques.

        Return
        -------
        list:
            - base (Origen X, Origen Y, Origen Z)
            - confBM (# Blocks X, # Blocks Y, # Blocks Z, Size Block, Height Block)
    """
    base = (50100.00, 46700.00, -225.0)
    confBM = (500, 190, 50, 10, 10)

class Rowindex(caractModel):
    """
    Devuelve la ubicación de un bloque en una tabla con dimensiones predefinidas.

        Args
        -----
            caractModel (List): Caracteristicas espaciales del modelo de bloques
        
        Parameters
        ----------
            point : list
                Type point (X, Y, Z)
        
        Returns
        -------
            int: 
                La posicion en la tabla
    """
    def __init__(self, point):
        self.point = point

    def idRow(self):
        numberRow = None

        def chCoor(base,size):
            return lambda coord:(((coord-base)-size/2)/size+1)

        i = chCoor(self.base[0],self.confBM[3])
        j = chCoor(self.base[1],self.confBM[3])
        k = chCoor(self.base[2],self.confBM[4])

        i_val = int(i(self.point[0]))
        j_val = int(j(self.point[1])-1)
        k_val = int(k(self.point[2])-1)
        numberRow = i_val + j_val*(self.confBM[0]) + k_val*(self.confBM[0]*self.confBM[1])

        return int(numberRow)
