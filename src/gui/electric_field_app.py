"""
Модуль реализует графический интерфейс для визуализации электрического поля.

Основной класс ElectricFieldApp создает окно с элементами управления и
областью для отображения силовых линий электрического поля. Позволяет
размещать заряды, создавать замкнутые области с распределенным зарядом
и визуализировать движение пробного заряда в поле.
"""

import numpy as np
import matplotlib.pyplot as plt
from tkinter import Tk, Button, Scale, HORIZONTAL, simpledialog
import matplotlib.patches as patches
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.path import Path

from ..models.body import Body
from ..models.charge import Charge
from ..models.polygon import Polygon


class ElectricFieldApp:
    """
    Основной класс приложения для визуализации электрического поля.
    
    Attributes:
        root (Tk): корневое окно приложения
        charges (list): список точечных зарядов
        polygons (list): список заряженных областей
        field_density (int): плотность силовых линий
        test_body (Body): пробное заряженное тело
    """
    def __init__(self, root):
        """
        Инициализация приложения.
        
        Args:
            root (Tk): корневое окно приложения
        """
        self.root = root
        self.root.title("Electric Field Visualization")

        self.charges = []
        self.polygons = []
        self.current_polygon = None
        self.field_density = 20
        self.test_body = None
        self.placing_test_charge = False
        self.animation_running = False

        # Создаем элементы управления
        self.start_button = Button(root, text="Старт визуализации", command=self.visualize_field)
        self.start_button.pack()

        self.create_polygon_button = Button(root, text="Создать замкнутую область", command=self.create_polygon)
        self.create_polygon_button.pack()

        self.place_test_charge_button = Button(root, text="Отметить пробный заряд",
                                               command=self.start_placing_test_charge)
        self.place_test_charge_button.pack()

        self.density_scale = Scale(root, from_=0.1, to=2.0, resolution=0.1, orient=HORIZONTAL, label="Частота линий")
        self.density_scale.set(1.0)
        self.density_scale.pack()

        # Создаем область для рисования
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        # Привязываем обработчики событий
        self.canvas.mpl_connect("button_press_event", self.add_point)
        self.canvas.mpl_connect("key_press_event", self.close_polygon)

        # Устанавливаем границы области
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)

        # Добавляем начальный заряд
        initial_charge = Charge(0.5, 0.5, 1)
        self.charges.append(initial_charge)
        self.ax.plot(0.5, 0.5, 'o', color='red')
        self.canvas.draw()

    def start_placing_test_charge(self):
        """Включает режим размещения пробного заряда."""
        self.placing_test_charge = True

    def add_point(self, event):
        """
        Обрабатывает добавление новой точки при клике мышью.
        
        Args:
            event: событие клика мыши
        """
        x, y = event.xdata, event.ydata
        if x is None or y is None:
            return

        if self.placing_test_charge:
            charge = simpledialog.askfloat("Заряд", "Введите величину заряда:", parent=self.root)
            if charge is not None:
                mass = simpledialog.askfloat("Масса", "Введите массу частицы:", parent=self.root)
                if mass is not None:
                    self.test_body = Body(x, y, mass, charge)
                    self.ax.plot(x, y, 'o', color='green')
                    self.canvas.draw()
                    self.placing_test_charge = False
                    self.animation_running = True
                    self.animate_test_charge()
        elif self.current_polygon:
            self.current_polygon.add_vertex(x, y)
            self.ax.plot(x, y, 'o', color='black')
            if self.current_polygon.completed:
                color = self.current_polygon.fill_color()
                polygon_patch = patches.Polygon(self.current_polygon.vertices, closed=True, color=color, alpha=0.3)
                self.ax.add_patch(polygon_patch)
                self.polygons.append(self.current_polygon)
                self.current_polygon = None
            self.canvas.draw()
        else:
            value = simpledialog.askfloat("Величина заряда", "Введите величину заряда:", parent=self.root)
            if value is not None:
                charge = Charge(x, y, value)
                self.charges.append(charge)
                color = 'red' if value > 0 else 'blue'
                self.ax.plot(x, y, 'o', color=color)
                self.canvas.draw()

    def create_polygon(self):
        """Создает новый многоугольник с заданными параметрами."""
        num_points = simpledialog.askinteger("Количество точек", "Введите количество точек области:", parent=self.root)
        charge_density = simpledialog.askfloat("Плотность заряда", "Введите плотность заряда:", parent=self.root)
        if num_points and charge_density:
            self.current_polygon = Polygon(charge_density, num_points)

    def close_polygon(self, event):
        """
        Завершает создание многоугольника при нажатии Enter.
        
        Args:
            event: событие нажатия клавиши
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
        Вычисляет суммарную напряженность поля от всех источников.
        
        Args:
            x (float): x-координата точки
            y (float): y-координата точки
            
        Returns:
            tuple: компоненты суммарного вектора напряженности (Ex, Ey)
        """
        Ex, Ey = 0, 0
        for charge in self.charges:
            ex, ey = charge.electric_field(x, y)
            Ex += ex
            Ey += ey
        for polygon in self.polygons:
            ex, ey = polygon.electric_field(x, y)
            Ex += ex
            Ey += ey
        return Ex, Ey

    def visualize_field(self):
        """Визуализирует силовые линии электрического поля."""
        # Настройка плотности сетки
        self.field_density = int(self.density_scale.get() * 20)
        
        # Создаем сетку точек
        x = np.linspace(0, 1, self.field_density)
        y = np.linspace(0, 1, self.field_density)
        X, Y = np.meshgrid(x, y)
        points = np.column_stack((X.ravel(), Y.ravel()))
        
        # Создаем маску для исключения точек внутри полигонов
        from matplotlib.path import Path
        mask = np.ones(len(points), dtype=bool)
        for polygon in self.polygons:
            path = Path(polygon.vertices)
            mask &= ~path.contains_points(points)
        
        # Вычисляем поле
        Ex = np.zeros(X.shape).ravel()
        Ey = np.zeros(Y.shape).ravel()
        valid_points = points[mask]
        if len(valid_points) > 0:
            field_x, field_y = zip(*[self.calculate_total_field(p[0], p[1]) for p in valid_points])
            Ex[mask] = field_x
            Ey[mask] = field_y
        
        Ex = Ex.reshape(X.shape)
        Ey = Ey.reshape(Y.shape)
        
        # Отрисовка
        self.ax.clear()
        E_magnitude = np.sqrt(Ex**2 + Ey**2)
        
        self.ax.streamplot(X, Y, Ex, Ey,
                          color=E_magnitude,
                          cmap='viridis',
                          density=1.2,
                          linewidth=1.2,
                          arrowsize=0.8,
                          arrowstyle='->',
                          broken_streamlines=False,
                          minlength=0.1)

        # Отрисовка зарядов и областей
        for charge in self.charges:
            color = 'red' if charge.value > 0 else 'blue'
            self.ax.plot(charge.x, charge.y, 'o', color=color, markersize=10)

        for polygon in self.polygons:
            polygon_patch = patches.Polygon(polygon.vertices, closed=True,
                                            color=polygon.fill_color(), alpha=0.3)
            self.ax.add_patch(polygon_patch)

        if self.test_body:
            self.ax.plot(self.test_body.x, self.test_body.y, 'o', color='green', markersize=8)

        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)
        self.canvas.draw()

    def animate_test_charge(self):
        """Анимирует движение пробного заряда в электрическом поле."""
        TIME_STEP = 0.001  # Шаг по времени для численного интегрирования
        MASS_SCALE = 100  # Коэффициент масштабирования массы
        FORCE_SCALE = 1e-8  # Коэффициент масштабирования силы
        
        # Рисуем начальное поле
        self.visualize_field()
        
        # Удаляем начальную точку пробного заряда
        if len(self.ax.lines) > 0:
            self.ax.lines[-1].remove()
        
        while self.animation_running:
            test_body_mass = self.test_body.mass/MASS_SCALE
            Ex, Ey = self.calculate_total_field(self.test_body.x, self.test_body.y)
            
            # Вычисляем силы
            Fx = Ex * self.test_body.charge * FORCE_SCALE
            Fy = Ey * self.test_body.charge * FORCE_SCALE

            self.test_body.update_position(Fx, Fy, TIME_STEP)

            # Обновляем отображение
            self.ax.plot(self.test_body.x, self.test_body.y, 'o', color='green', markersize=8)
            self.canvas.draw()
            
            # Удаляем старую точку
            self.ax.lines[-1].remove()
            
            # Проверяем выход за границы
            if (self.test_body.x < 0 or self.test_body.x > 1 or 
                self.test_body.y < 0 or self.test_body.y > 1):
                self.animation_running = False
                return
                
            self.root.update()