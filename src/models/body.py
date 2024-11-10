class Body:
    """Класс, представляющий физическое тело с массой и зарядом."""
    
    def __init__(self, x: float, y: float, mass: float, charge: float):
        self.x = x
        self.y = y
        self.mass = mass
        self.charge = charge
        self.velocity_x = 0
        self.velocity_y = 0

    def update_position(self, force_x: float, force_y: float, time_step: float):
        acceleration_x = force_x / self.mass
        acceleration_y = force_y / self.mass
        
        self.velocity_x += acceleration_x * time_step
        self.velocity_y += acceleration_y * time_step
        
        self.x += self.velocity_x * time_step
        self.y += self.velocity_y * time_step 