import numpy as np
import matplotlib.pyplot as plt
from tkinter import Tk, Button, Scale, HORIZONTAL, simpledialog
import matplotlib.patches as patches
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


class Body:
    """
    Класс, представляющий физическое тело с массой и зарядом.
    
    Attributes:
        x (float): X-координата тела
        y (float): Y-координата тела 
        mass (float): Масса тела
        charge (float): Заряд тела
        velocity_x (float): Скорость тела по оси X
        velocity_y (float): Скорость тела по оси Y
    """
    def __init__(self, x: float, y: float, mass: float, charge: float):
        """
        Инициализация тела.

        Args:
            x (float): Начальная X-координата
            y (float): Начальная Y-координата
            mass (float): Масса тела
            charge (float): Заряд тела
        """
        self.x = x
        self.y = y
        self.mass = mass
        self.charge = charge
        self.velocity_x = 0
        self.velocity_y = 0

    def update_position(self, force_x: float, force_y: float, time_step: float):
        """
        Обновляет положение и скорость тела под действием сил.

        Args:
            force_x (float): Сила по оси X
            force_y (float): Сила по оси Y
            time_step (float): Шаг времени для расчета
        """
        # Расчет ускорения по второму закону Ньютона: a = F/m
        acceleration_x = force_x / self.mass
        acceleration_y = force_y / self.mass

        # Обновление скорости: v = v0 + at
        self.velocity_x += acceleration_x * time_step
        self.velocity_y += acceleration_y * time_step

        # Обновление положения: x = x0 + vt
        self.x += self.velocity_x * time_step
        self.y += self.velocity_y * time_step


class Charge:
    """
    Класс, представляющий точечный электрический заряд.
    
    Attributes:
        x (float): X-координата заряда
        y (float): Y-координата заряда
        value (float): Величина заряда в Кулонах
    """
    def __init__(self, x: float, y: float, value: float):
        """
        Инициализация заряда.

        Args:
            x (float): X-координата заряда
            y (float): Y-координата заряда
            value (float): Величина заряда в Кулонах
        """
        self.x = x
        self.y = y
        self.value = value

    def electric_field(self, point_x: float, point_y: float) -> tuple[float, float]:
        """
        Вычисляет напряженность электрического поля в заданной точке.

        Args:
            point_x (float): X-координата точки
            point_y (float): Y-координата точки

        Returns:
            tuple[float, float]: Компоненты вектора напряженности (Ex, Ey)
        """
        COULOMB_CONSTANT = 8.99e9  # Постоянная Кулона в Н·м²/Кл²
        MIN_DISTANCE = 0.001  # Минимальное расстояние для избежания деления на ноль
        
        # Вычисление расстояния от заряда до точки
        distance_x = point_x - self.x
        distance_y = point_y - self.y
        distance = np.sqrt(distance_x ** 2 + distance_y ** 2)
        
        if distance < MIN_DISTANCE:
            return 0, 0
            
        # Расчет напряженности по закону Кулона
        field_magnitude = COULOMB_CONSTANT * self.value / distance ** 2
        field_x = field_magnitude * distance_x / distance
        field_y = field_magnitude * distance_y / distance
        
        return field_x, field_y


class Polygon:
    """
    Класс, представляющий многоугольную область с распределенным зарядом.
    
    Attributes:
        charge_density (float): Плотность заряда на единицу площади
        num_points (int): Количество вершин многоугольника
        vertices (list): Список координат вершин
        completed (bool): Флаг завершенности построения многоугольника
    """
    def __init__(self, charge_density: float, num_points: int):
        """
        Инициализация многоугольника.

        Args:
            charge_density (float): Плотность заряда
            num_points (int): Количество вершин
        """
        self.charge_density = charge_density
        self.num_points = num_points
        self.vertices = []
        self.completed = False

    def add_vertex(self, x: float, y: float):
        """
        Добавляет новую вершину в многоугольник.

        Args:
            x (float): X-координата вершины
            y (float): Y-координата вершины
        """
        self.vertices.append((x, y))
        if len(self.vertices) == self.num_points:
            self.completed = True

    def electric_field(self, point_x: float, point_y: float) -> tuple[float, float]:
        """
        Вычисляет напряженность электрического поля от многоугольника в точке.

        Args:
            point_x (float): X-координата точки
            point_y (float): Y-координата точки

        Returns:
            tuple[float, float]: Компоненты вектора напряженности (Ex, Ey)
        """
        field_x, field_y = 0, 0
        COULOMB_CONSTANT = 8.99e9
        GRID_SIZE = 10  # Количество разбиений по каждой оси
        MIN_DISTANCE = 0.001  # Минимальное расстояние для расчетов
        
        # Находим границы полигона
        x_min = min(vertex[0] for vertex in self.vertices)
        x_max = max(vertex[0] for vertex in self.vertices)
        y_min = min(vertex[1] for vertex in self.vertices)
        y_max = max(vertex[1] for vertex in self.vertices)
        
        # Размеры элементарной ячейки
        cell_width = (x_max - x_min) / GRID_SIZE
        cell_height = (y_max - y_min) / GRID_SIZE
        cell_area = cell_width * cell_height
        
        from matplotlib.path import Path
        polygon_path = Path(self.vertices)
        
        # Расчет поля от каждой ячейки
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
        """
        Определяет цвет заливки многоугольника в зависимости от знака заряда.

        Returns:
            str: Цвет заливки ('red' для положительного заряда, 'blue' для отрицательного)
        """
        return 'red' if self.charge_density > 0 else 'blue'


class ElectricFieldApp:
    """
    Класс для визуализации электрического поля и взаимодействия с пользователем.
    
    Attributes:
        root (Tk): Корневое окно приложения
        charges (list): Список точечных зарядов
        polygons (list): Список полигонов с распределенным зарядом
        current_polygon (Polygon): Текущий создаваемый полигон
        field_density (int): Плотность сетки для визуализации поля
        test_body (Body): Пробное тело для анимации движения
        placing_test_charge (bool): Флаг размещения пробного заряда
        animation_running (bool): Флаг работы анимации
    """
    def __init__(self, root):
        """
        Инициализация приложения.
        
        Args:
            root (Tk): Корневое окно приложения
        """
        self.root = root
        self.root.title("Визуализация электрического поля")

        # Инициализация основных переменных
        self.charges = []
        self.polygons = []
        self.current_polygon = None
        self.field_density = 20  # Базовая плотность сетки
        self.test_body = None
        self.placing_test_charge = False
        self.animation_running = False

        # Создание элементов управления
        self._create_control_elements()
        
        # Настройка matplotlib
        self._setup_matplotlib()
        
        # Добавление начального заряда
        self._add_initial_charge()

    def _create_control_elements(self):
        """Создание кнопок и элементов управления"""
        self.start_button = Button(self.root, text="Старт визуализации", command=self.visualize_field)
        self.start_button.pack()

        self.create_polygon_button = Button(self.root, text="Создать замкнутую область", command=self.create_polygon)
        self.create_polygon_button.pack()

        self.place_test_charge_button = Button(self.root, text="Отметить пробный заряд",
                                             command=self.start_placing_test_charge)
        self.place_test_charge_button.pack()

        self.density_scale = Scale(self.root, from_=0.1, to=2.0, resolution=0.1, 
                                 orient=HORIZONTAL, label="Частота линий")
        self.density_scale.set(1.0)
        self.density_scale.pack()

    def _setup_matplotlib(self):
        """Настройка области рисования matplotlib"""
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        # Привязка обработчиков событий
        self.canvas.mpl_connect("button_press_event", self.add_point)
        self.canvas.mpl_connect("key_press_event", self.close_polygon)

        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)

    def _add_initial_charge(self):
        """Добавление начального положительного заряда в центр"""
        initial_charge = Charge(0.5, 0.5, 1)  # Заряд величиной 1 в центре
        self.charges.append(initial_charge)
        self.ax.plot(0.5, 0.5, 'o', color='red')
        self.canvas.draw()

    def start_placing_test_charge(self):
        """Включение режима размещения пробного заряда"""
        self.placing_test_charge = True

    def add_point(self, event):
        """
        Обработка клика мыши для добавления объектов.
        
        Args:
            event: Событие клика мыши
        """
        click_x, click_y = event.xdata, event.ydata
        if click_x is None or click_y is None:
            return

        if self.placing_test_charge:
            self._handle_test_charge_placement(click_x, click_y)
        elif self.current_polygon:
            self._handle_polygon_point(click_x, click_y)
        else:
            self._handle_charge_placement(click_x, click_y)

    def _handle_test_charge_placement(self, x, y):
        """Обработка размещения пробного заряда"""
        charge_value = simpledialog.askfloat("Заряд", "Введите величину заряда:", parent=self.root)
        if charge_value is not None:
            mass_value = simpledialog.askfloat("Масса", "Введите массу частицы:", parent=self.root)
            if mass_value is not None:
                self.test_body = Body(x, y, mass_value, charge_value)
                self.ax.plot(x, y, 'o', color='green')
                self.canvas.draw()
                self.placing_test_charge = False
                self.animation_running = True
                self.animate_test_charge()

    def _handle_polygon_point(self, x, y):
        """Обработка добавления точки в полигон"""
        self.current_polygon.add_vertex(x, y)
        self.ax.plot(x, y, 'o', color='black')
        if self.current_polygon.completed:
            color = self.current_polygon.fill_color()
            polygon_patch = patches.Polygon(self.current_polygon.vertices, closed=True, color=color, alpha=0.3)
            self.ax.add_patch(polygon_patch)
            self.polygons.append(self.current_polygon)
            self.current_polygon = None
        self.canvas.draw()

    def _handle_charge_placement(self, x, y):
        """Обработка размещения точечного заряда"""
        charge_value = simpledialog.askfloat("Величина заряда", "Введите величину заряда:", parent=self.root)
        if charge_value is not None:
            charge = Charge(x, y, charge_value)
            self.charges.append(charge)
            color = 'red' if charge_value > 0 else 'blue'
            self.ax.plot(x, y, 'o', color=color)
            self.canvas.draw()

    def create_polygon(self):
        """Создание нового полигона"""
        num_points = simpledialog.askinteger("Количество точек", "Введите количество точек области:", parent=self.root)
        charge_density = simpledialog.askfloat("Плотность заряда", "Введите плотность заряда:", parent=self.root)
        if num_points and charge_density:
            self.current_polygon = Polygon(charge_density, num_points)

    def close_polygon(self, event):
        """
        Завершение создания полигона по нажатию Enter.
        
        Args:
            event: Событие нажатия клавиши
        """
        if event.key == 'enter' and self.current_polygon and self.current_polygon.completed:
            color = self.current_polygon.fill_color()
            polygon_patch = patches.Polygon(self.current_polygon.vertices, closed=True, color=color, alpha=0.3)
            self.ax.add_patch(polygon_patch)
            self.polygons.append(self.current_polygon)
            self.current_polygon = None
            self.canvas.draw()

    def calculate_total_field(self, x, y):
        """
        Расчет суммарного электрического поля в точке.
        
        Args:
            x (float): X-координата точки
            y (float): Y-координата точки
            
        Returns:
            tuple: Компоненты вектора напряженности (Ex, Ey)
        """
        field_x, field_y = 0, 0
        for charge in self.charges:
            ex, ey = charge.electric_field(x, y)
            field_x += ex
            field_y += ey
        for polygon in self.polygons:
            ex, ey = polygon.electric_field(x, y)
            field_x += ex
            field_y += ey
        return field_x, field_y

    def visualize_field(self):
        """Визуализация силовых линий электрического поля"""
        # Настройка плотности сетки
        self.field_density = int(self.density_scale.get() * 20)
        
        # Создание сетки точек
        x_coords = np.linspace(0, 1, self.field_density)
        y_coords = np.linspace(0, 1, self.field_density)
        X, Y = np.meshgrid(x_coords, y_coords)
        points = np.column_stack((X.ravel(), Y.ravel()))
        
        # Создание маски для точек вне полигонов
        from matplotlib.path import Path
        mask = np.ones(len(points), dtype=bool)
        for polygon in self.polygons:
            path = Path(polygon.vertices)
            mask &= ~path.contains_points(points)
        
        # Вычисление поля
        Ex = np.zeros(X.shape).ravel()
        Ey = np.zeros(Y.shape).ravel()
        valid_points = points[mask]
        if len(valid_points) > 0:
            field_x, field_y = zip(*[self.calculate_total_field(p[0], p[1]) for p in valid_points])
            Ex[mask] = field_x
            Ey[mask] = field_y
        
        Ex = Ex.reshape(X.shape)
        Ey = Ey.reshape(Y.shape)
        
        self._draw_field(X, Y, Ex, Ey)

    def _draw_field(self, X, Y, Ex, Ey):
        """Отрисовка силовых линий и объектов"""
        self.ax.clear()
        field_magnitude = np.sqrt(Ex**2 + Ey**2)
        
        # Отрисовка силовых линий
        self.ax.streamplot(X, Y, Ex, Ey,
                          color=field_magnitude,
                          cmap='viridis',
                          density=1,
                          linewidth=1.2,
                          arrowsize=0.8,
                          arrowstyle='->',
                          broken_streamlines=False,
                          minlength=0.1)

        # Отрисовка зарядов
        for charge in self.charges:
            color = 'red' if charge.value > 0 else 'blue'
            self.ax.plot(charge.x, charge.y, 'o', color=color, markersize=10)

        # Отрисовка полигонов
        for polygon in self.polygons:
            polygon_patch = patches.Polygon(polygon.vertices, closed=True,
                                         color=polygon.fill_color(), alpha=0.3)
            self.ax.add_patch(polygon_patch)

        # Отрисовка пробного заряда
        if self.test_body:
            self.ax.plot(self.test_body.x, self.test_body.y, 'o', color='green', markersize=8)

        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)
        self.canvas.draw()

    def animate_test_charge(self):
        """Анимация движения пробного заряда"""
        TIME_STEP = 0.001  # Шаг по времени для расчета движения
        FORCE_SCALE = 1e-8  # Масштабный коэффициент для сил
        
        # Начальная визуализация поля
        self.visualize_field()
        
        # Удаление начальной точки пробного заряда
        if len(self.ax.lines) > 0:
            self.ax.lines[-1].remove()
        
        while self.animation_running:
            # Расчет массы и сил
            scaled_mass = self.test_body.mass/100
            field_x, field_y = self.calculate_total_field(self.test_body.x, self.test_body.y)
            force_x = field_x * self.test_body.charge * FORCE_SCALE
            force_y = field_y * self.test_body.charge * FORCE_SCALE

            # Обновление положения
            self.test_body.update_position(force_x, force_y, TIME_STEP)

            # Отрисовка текущего положения
            self.ax.plot(self.test_body.x, self.test_body.y, 'o', color='green', markersize=8)
            self.canvas.draw()
            
            # Очистка предыдущего положения
            self.ax.lines[-1].remove()
            
            # Проверка выхода за границы
            if not (0 <= self.test_body.x <= 1 and 0 <= self.test_body.y <= 1):
                self.animation_running = False
                return
                
            self.root.update()


if __name__ == "__main__":
    root = Tk()
    app = ElectricFieldApp(root)
    root.mainloop()
