import numpy as np
from matplotlib.path import Path
from typing import List, Tuple, Optional

class Polygon:
    """
    Класс для представления заряженной области в форме многоугольника.
    
    Attributes:
        charge_density (float): плотность заряда
        num_points (int): количество вершин многоугольника
        vertices (list): список координат вершин
        completed (bool): завершен ли многоугольник
    """
    def __init__(self, charge_density: float, num_points: int) -> None:
        """
        Инициализация многоугольника.
        
        Args:
            charge_density (float): плотность заряда
            num_points (int): количество вершин
        """
        self.charge_density: float = charge_density
        self.num_points: int = num_points
        self.vertices: List[Tuple[float, float]] = []
        self.completed: bool = False

    def add_vertex(self, x: float, y: float) -> None:
        """
        Добавляет новую вершину в многоугольник.
        
        Args:
            x (float): x-координата вершины
            y (float): y-координата вершины
        """
        self.vertices.append((x, y))
        if len(self.vertices) == self.num_points:
            self.completed = True

    def electric_field(self, x: float, y: float) -> Tuple[float, float]:
        """
        Вычисляет напряженность электрического поля от заряженной области.
        
        Args:
            x (float): x-координата точки
            y (float): y-координата точки
            
        Returns:
            tuple: компоненты вектора напряженности (Ex, Ey)
        """
        Ex, Ey = 0.0, 0.0
        COULOMB_CONSTANT: float = 8.99e9
        GRID_SIZE: int = 10  # Количество разбиений по каждой оси
        MIN_DISTANCE: float = 0.001  # Минимальное расстояние для расчетов
        
        # Находим границы полигона
        x_min: float = min(v[0] for v in self.vertices)
        x_max: float = max(v[0] for v in self.vertices)
        y_min: float = min(v[1] for v in self.vertices)
        y_max: float = max(v[1] for v in self.vertices)
        
        # Размеры элементарной ячейки
        dx: float = (x_max - x_min) / GRID_SIZE
        dy: float = (y_max - y_min) / GRID_SIZE
        dA: float = dx * dy  # Площадь элементарной ячейки
        
        from matplotlib.path import Path
        polygon_path: Path = Path(self.vertices)
        
        # Суммируем вклады от всех элементарных ячеек
        for i in range(GRID_SIZE):
            for j in range(GRID_SIZE):
                xi: float = x_min + (i + 0.5) * dx
                yi: float = y_min + (j + 0.5) * dy
                
                if polygon_path.contains_point((xi, yi)):
                    r_x: float = x - xi
                    r_y: float = y - yi
                    r: float = np.sqrt(r_x**2 + r_y**2)
                    
                    if r > MIN_DISTANCE:
                        dq: float = self.charge_density * dA
                        e: float = COULOMB_CONSTANT * dq / r**2
                        Ex += e * r_x / r
                        Ey += e * r_y / r
        
        return Ex, Ey

    def fill_color(self) -> str:
        """
        Определяет цвет заливки многоугольника в зависимости от знака заряда.
        
        Returns:
            str: цвет заливки ('red' для положительного заряда, 'blue' для отрицательного)
        """
        return 'red' if self.charge_density > 0 else 'blue'