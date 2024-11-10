import numpy as np

class Charge:
    """
    Класс для представления точечного заряда.
    
    Attributes:
        x (float): x-координата заряда
        y (float): y-координата заряда
        value (float): величина заряда в Кулонах
    """
    def __init__(self, x, y, value):
        """
        Инициализация точечного заряда.
        
        Args:
            x (float): x-координата
            y (float): y-координата
            value (float): величина заряда
        """
        self.x = x
        self.y = y
        self.value = value

    def electric_field(self, x, y):
        """
        Вычисляет напряженность электрического поля в заданной точке.
        
        Args:
            x (float): x-координата точки
            y (float): y-координата точки
            
        Returns:
            tuple: компоненты вектора напряженности (Ex, Ey)
        """
        COULOMB_CONSTANT = 8.99e9  # Постоянная Кулона в Н·м²/Кл²
        MIN_DISTANCE = 0.001  # Минимальное расстояние для избежания деления на ноль
        
        r_x = x - self.x
        r_y = y - self.y
        r = np.sqrt(r_x ** 2 + r_y ** 2)
        
        if r < MIN_DISTANCE:
            return 0, 0
            
        e = COULOMB_CONSTANT * self.value / r ** 2
        ex = e * r_x / r
        ey = e * r_y / r
        return ex, ey