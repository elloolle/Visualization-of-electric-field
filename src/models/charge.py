import numpy as np

class Charge:
    """Класс, представляющий точечный электрический заряд."""
    
    def __init__(self, x: float, y: float, value: float):
        self.x = x
        self.y = y
        self.value = value

    def electric_field(self, point_x: float, point_y: float) -> tuple[float, float]:
        COULOMB_CONSTANT = 8.99e9
        MIN_DISTANCE = 0.001
        
        distance_x = point_x - self.x
        distance_y = point_y - self.y
        distance = np.sqrt(distance_x ** 2 + distance_y ** 2)
        
        if distance < MIN_DISTANCE:
            return 0, 0
            
        field_magnitude = COULOMB_CONSTANT * self.value / distance ** 2
        field_x = field_magnitude * distance_x / distance
        field_y = field_magnitude * distance_y / distance
        
        return field_x, field_y 