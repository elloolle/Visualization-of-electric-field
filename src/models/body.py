class Body:
    """
    Класс для представления движущегося заряженного тела.
    
    Attributes:
        x (float): x-координата тела
        y (float): y-координата тела 
        mass (float): масса тела в кг
        charge (float): заряд тела в Кулонах
        vx (float): скорость по x в м/с
        vy (float): скорость по y в м/с
    """
    def __init__(self, x: float, y: float, mass: float, charge: float) -> None:
        """
        Инициализация тела.
        
        Args:
            x (float): начальная x-координата
            y (float): начальная y-координата
            mass (float): масса тела
            charge (float): заряд тела
        """
        self.x: float = x
        self.y: float = y
        self.mass: float = mass
        self.charge: float = charge
        self.vx: float = 0  # Начальная скорость по x
        self.vy: float = 0  # Начальная скорость по y

    def update_position(self, fx: float, fy: float, dt: float) -> None:
        """
        Обновляет положение и скорость тела согласно законам Ньютона.
        
        Args:
            fx (float): сила по x
            fy (float): сила по y
            dt (float): временной шаг
        """
        # Вычисляем ускорение по второму закону Ньютона: F = ma
        ax: float = fx / self.mass
        ay: float = fy / self.mass

        # Обновляем скорости
        self.vx += ax * dt
        self.vy += ay * dt

        # Обновляем координаты
        self.x += self.vx * dt
        self.y += self.vy * dt