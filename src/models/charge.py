import numpy as np
from typing import Tuple

class Charge:
    """
    Класс для представления точечного заряда.
    
    Attributes:
        x (float): x-координата заряда
        y (float): y-координата заряда
        value (float): величина заряда в Кулонах
    """
    def __init__(self, x: float, y: float, value: float) -> None:
        """
        Инициализация точечного заряда.
        
        Args:
            x (float): x-координата
            y (float): y-координата
            value (float): величина заряда
        """
        self.x: float = x
        self.y: float = y
        self.value: float = value

    def electric_field(self, x: float, y: float) -> Tuple[float, float]:
        """
        Вычисляет напряженность электрического поля в заданной точке.
        
        Args:
            x (float): x-координата точки
            y (float): y-координата точки
            
        Returns:
            tuple: компоненты вектора напряженности (Ex, Ey)
        """
        COULOMB_CONSTANT: float = 8.99e9  # Постоянная Кулона в Н·м²/Кл²
        MIN_DISTANCE: float = 0.001  # Минимальное расстояние для избежания деления на ноль
        
        r_x: float = x - self.x
        r_y: float = y - self.y
        r: float = np.sqrt(r_x ** 2 + r_y ** 2)
        
        if r < MIN_DISTANCE:
            return 0, 0
            
        e: float = COULOMB_CONSTANT * self.value / r ** 2
        ex: float = e * r_x / r
        ey: float = e * r_y / r
        return ex, ey