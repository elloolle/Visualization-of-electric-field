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
    def __init__(self, x, y, mass, charge):
        """
        Инициализация тела.
        
        Args:
            x (float): начальная x-координата
            y (float): начальная y-координата
            mass (float): масса тела
            charge (float): заряд тела
        """
        self.x = x
        self.y = y
        self.mass = mass
        self.charge = charge
        self.vx = 0  # Начальная скорость по x
        self.vy = 0  # Начальная скорость по y

    def update_position(self, fx, fy, dt):
        """
        Обновляет положение и скорость тела согласно законам Ньютона.
        
        Args:
            fx (float): сила по x
            fy (float): сила по y
            dt (float): временной шаг
        """
        # Вычисляем ускорение по второму закону Ньютона: F = ma
        ax = fx / self.mass
        ay = fy / self.mass

        # Обновляем скорости
        self.vx += ax * dt
        self.vy += ay * dt

        # Обновляем координаты
        self.x += self.vx * dt
        self.y += self.vy * dt