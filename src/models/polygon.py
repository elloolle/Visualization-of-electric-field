import numpy as np
from matplotlib.path import Path

class Polygon:
    """
    Класс для представления заряженной области в форме многоугольника.
    
    Attributes:
        charge_density (float): плотность заряда
        num_points (int): количество вершин многоугольника
        vertices (list): список координат вершин
        completed (bool): завершен ли многоугольник
    """
    def __init__(self, charge_density, num_points):
        """
        Инициализация многоугольника.
        
        Args:
            charge_density (float): плотность заряда
            num_points (int): количество вершин
        """
        self.charge_density = charge_density
        self.num_points = num_points
        self.vertices = []
        self.completed = False

    def add_vertex(self, x, y):
        """
        Добавляет новую вершину в многоугольник.
        
        Args:
            x (float): x-координата вершины
            y (float): y-координата вершины
        """
        self.vertices.append((x, y))
        if len(self.vertices) == self.num_points:
            self.completed = True

    def electric_field(self, x, y):
        """
        Вычисляет напряженность электрического поля от заряженной области.
        
        Args:
            x (float): x-координата точки
            y (float): y-координата точки
            
        Returns:
            tuple: компоненты вектора напряженности (Ex, Ey)
        """
        Ex, Ey = 0, 0
        COULOMB_CONSTANT = 8.99e9
        GRID_SIZE = 10  # Количество разбиений по каждой оси
        MIN_DISTANCE = 0.001  # Минимальное расстояние для расчетов
        
        # Находим границы полигона
        x_min = min(v[0] for v in self.vertices)
        x_max = max(v[0] for v in self.vertices)
        y_min = min(v[1] for v in self.vertices)
        y_max = max(v[1] for v in self.vertices)
        
        # Размеры элементарной ячейки
        dx = (x_max - x_min) / GRID_SIZE
        dy = (y_max - y_min) / GRID_SIZE
        dA = dx * dy  # Площадь элементарной ячейки
        
        from matplotlib.path import Path
        polygon_path = Path(self.vertices)
        
        # Суммируем вклады от всех элементарных ячеек
        for i in range(GRID_SIZE):
            for j in range(GRID_SIZE):
                xi = x_min + (i + 0.5) * dx
                yi = y_min + (j + 0.5) * dy
                
                if polygon_path.contains_point((xi, yi)):
                    r_x = x - xi
                    r_y = y - yi
                    r = np.sqrt(r_x**2 + r_y**2)
                    
                    if r > MIN_DISTANCE:
                        dq = self.charge_density * dA
                        e = COULOMB_CONSTANT * dq / r**2
                        Ex += e * r_x / r
                        Ey += e * r_y / r
        
        return Ex, Ey

    def fill_color(self):
        """
        Определяет цвет заливки многоугольника в зависимости от знака заряда.
        
        Returns:
            str: цвет заливки ('red' для положительного заряда, 'blue' для отрицательного)
        """
        return 'red' if self.charge_density > 0 else 'blue'