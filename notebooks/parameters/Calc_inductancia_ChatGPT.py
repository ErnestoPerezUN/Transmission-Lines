#Códico para calcular la inductancia y capacitancia de una línea de transmisión
# CODIGO REALIZADO POR CHATGPT
# prompt: "Crear un código en python para calcular 
#           los parámetros de inductancia y capacitancia 
#           de múltiples conductores aéreos, teniendo 
#           una estructura de objetos"
# fecha: 2025-05-22
# Comentarios:
#   El código aparentemente no tiene errores, pero no se ha probado


import math

# Constantes
MU_0 = 4 * math.pi * 1e-7       # Permeabilidad del vacío (H/m)
EPSILON_0 = 8.854e-12           # Permitividad del vacío (F/m)

class Conductor:
    def __init__(self, name, x, y, radius):
        self.name = name
        self.x = x
        self.y = y
        self.radius = radius
    
    def distance_to(self, other):
        return math.hypot(self.x - other.x, self.y - other.y)

class TransmissionLine:
    def __init__(self, conductors):
        self.conductors = conductors

    def calc_GMD(self):
        """Distancia media geométrica (GMD) entre conductores"""
        n = len(self.conductors)
        product = 1.0
        count = 0

        for i in range(n):
            for j in range(i+1, n):
                d = self.conductors[i].distance_to(self.conductors[j])
                product *= d
                count += 1

        return product ** (1 / count)

    def calc_GMR(self):
        """Radio medio geométrico (GMR) de todos los conductores"""
        product = 1.0
        for conductor in self.conductors:
            # r' = 0.7788 * r para un solo hilo
            product *= 0.7788 * conductor.radius
        return product ** (1 / len(self.conductors))

    def inductance_per_meter(self):
        GMD = self.calc_GMD()
        GMR = self.calc_GMR()
        return 2e-7 * math.log(GMD / GMR)  # H/m

    def capacitance_per_meter(self):
        GMD = self.calc_GMD()
        r_mean = sum(c.radius for c in self.conductors) / len(self.conductors)
        return (2 * math.pi * EPSILON_0) / math.log(GMD / r_mean)  # F/m

    def summary(self):
        L = self.inductance_per_meter()
        C = self.capacitance_per_meter()
        print("Inductancia total (H/m):", L)
        print("Capacitancia total (F/m):", C)

# Ejemplo de uso
if __name__ == "__main__":
    conductores = [
        Conductor("A", 0, 0, 0.01),
        Conductor("B", 0.3, 0, 0.01),
        Conductor("C", 0.15, 0.26, 0.01)
    ]

    linea = TransmissionLine(conductores)
    linea.summary()
