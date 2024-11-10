import numpy as np
from matplotlib.path import Path

class Polygon:
    """Класс, представляющий многоугольную область с распределенным зарядом."""
    
    def __init__(self, charge_density: float, num_points: int):
        self.charge_density = charge_density
        self.num_points = num_points
        self.vertices = []
        self.completed = False

    def add_vertex(self, x: float, y: float):
        self.vertices.append((x, y))
        if len(self.vertices) == self.num_points:
            self.completed = True

    def electric_field(self, point_x: float, point_y: float) -> tuple[float, float]:
        field_x, field_y = 0, 0
        COULOMB_CONSTANT = 8.99e9
        GRID_SIZE = 10
        MIN_DISTANCE = 0.001
        
        x_min = min(vertex[0] for vertex in self.vertices)
        x_max = max(vertex[0] for vertex in self.vertices)
        y_min = min(vertex[1] for vertex in self.vertices)
        y_max = max(vertex[1] for vertex in self.vertices)
        
        cell_width = (x_max - x_min) / GRID_SIZE
        cell_height = (y_max - y_min) / GRID_SIZE
        cell_area = cell_width * cell_height
        
        polygon_path = Path(self.vertices)
        
        for i in range(GRID_SIZE):
            for j in range(GRID_SIZE):
                cell_center_x = x_min + (i + 0.5) * cell_width
                cell_center_y = y_min + (j + 0.5) * cell_height
                
                if polygon_path.contains_point((cell_center_x, cell_center_y)):
                    distance_x = point_x - cell_center_x
                    distance_y = point_y - cell_center_y
                    distance = np.sqrt(distance_x**2 + distance_y**2)
                    
                    if distance > MIN_DISTANCE:
                        cell_charge = self.charge_density * cell_area
                        field_magnitude = COULOMB_CONSTANT * cell_charge / distance**2
                        field_x += field_magnitude * distance_x / distance
                        field_y += field_magnitude * distance_y / distance
        
        return field_x, field_y

    def fill_color(self) -> str:
        return 'red' if self.charge_density > 0 else 'blue' 