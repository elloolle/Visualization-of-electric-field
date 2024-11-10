import numpy as np
import matplotlib.pyplot as plt
from tkinter import Tk, Button, Scale, HORIZONTAL, simpledialog
import matplotlib.patches as patches
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import time


class Body:
    def __init__(self, x, y, mass, charge):
        self.x = x
        self.y = y
        self.mass = mass
        self.charge = charge
        self.vx = 0
        self.vy = 0

    def update_position(self, fx, fy, dt):
        # Обновляем скорость и положение по законам Ньютона
        ax = fx / self.mass
        ay = fy / self.mass

        self.vx += ax * dt
        self.vy += ay * dt

        self.x += self.vx * dt
        self.y += self.vy * dt


class Charge:
    def __init__(self, x, y, value):
        self.x = x
        self.y = y
        self.value = value

    def electric_field(self, x, y):
        k = 8.99e9  # постоянная электростатического взаимодействия (в Н·м²/Кл²)
        r_x = x - self.x
        r_y = y - self.y
        r = np.sqrt(r_x ** 2 + r_y ** 2)
        if r < 0.001:  # Избегаем деления на ноль
            return 0, 0
        e = k * self.value / r ** 2
        ex = e * r_x / r
        ey = e * r_y / r
        return ex, ey


class Polygon:
    def __init__(self, charge_density, num_points):
        self.charge_density = charge_density
        self.num_points = num_points
        self.vertices = []
        self.completed = False

    def add_vertex(self, x, y):
        self.vertices.append((x, y))
        if len(self.vertices) == self.num_points:
            self.completed = True

    def electric_field(self, x, y):
        Ex, Ey = 0, 0
        k = 8.99e9
        
        # Находим границы полигона
        x_min = min(v[0] for v in self.vertices)
        x_max = max(v[0] for v in self.vertices)
        y_min = min(v[1] for v in self.vertices)
        y_max = max(v[1] for v in self.vertices)
        
        # Разбиваем на сетку
        grid_size = 10  # количество разбиений по каждой оси
        dx = (x_max - x_min) / grid_size
        dy = (y_max - y_min) / grid_size
        
        # Площадь элементарной ячейки
        dA = dx * dy
        
        from matplotlib.path import Path
        polygon_path = Path(self.vertices)
        
        # Перебираем все ячейки сетки
        for i in range(grid_size):
            for j in range(grid_size):
                # Координаты центра текущей ячейки
                xi = x_min + (i + 0.5) * dx
                yi = y_min + (j + 0.5) * dy
                
                # Проверяем, находится ли точка внутри полигона
                if polygon_path.contains_point((xi, yi)):
                    # Расчет поля от этой ячейки
                    r_x = x - xi
                    r_y = y - yi
                    r = np.sqrt(r_x**2 + r_y**2)
                    
                    if r > 0.001:  # Избегаем деления на ноль
                        # Заряд ячейки = плотность * площадь
                        dq = self.charge_density * dA
                        e = k * dq / r**2
                        Ex += e * r_x / r
                        Ey += e * r_y / r
        
        return Ex, Ey

    def fill_color(self):
        return 'red' if self.charge_density > 0 else 'blue'


class ElectricFieldApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Electric Field Visualization")

        self.charges = []
        self.polygons = []
        self.current_polygon = None
        self.field_density = 20
        self.test_body = None
        self.placing_test_charge = False
        self.animation_running = False

        # Элементы управления
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

        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        self.canvas.mpl_connect("button_press_event", self.add_point)
        self.canvas.mpl_connect("key_press_event", self.close_polygon)

        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)

        initial_charge = Charge(0.5, 0.5, 1)
        self.charges.append(initial_charge)
        self.ax.plot(0.5, 0.5, 'o', color='red')
        self.canvas.draw()

    def start_placing_test_charge(self):
        self.placing_test_charge = True

    def add_point(self, event):
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
        num_points = simpledialog.askinteger("Количество точек", "Введите количество точек области:", parent=self.root)
        charge_density = simpledialog.askfloat("Плотность заряда", "Введите плотность заряда:", parent=self.root)
        if num_points and charge_density:
            self.current_polygon = Polygon(charge_density, num_points)

    def close_polygon(self, event):
        if event.key == 'enter' and self.current_polygon and self.current_polygon.completed:
            color = self.current_polygon.fill_color()
            polygon_patch = patches.Polygon(self.current_polygon.vertices, closed=True, color=color, alpha=0.3)
            self.ax.add_patch(polygon_patch)
            self.polygons.append(self.current_polygon)
            self.current_polygon = None
            self.canvas.draw()

    def calculate_total_field(self, x, y):
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
        # Уменьшаем плотность сетки
        self.field_density = int(self.density_scale.get() * 20)  # уменьшено с 30
        
        # Используем векторизацию numpy вместо циклов
        x = np.linspace(0, 1, self.field_density)
        y = np.linspace(0, 1, self.field_density)
        X, Y = np.meshgrid(x, y)
        points = np.column_stack((X.ravel(), Y.ravel()))
        
        # Создаем маску для точек внутри полигонов
        from matplotlib.path import Path
        mask = np.ones(len(points), dtype=bool)
        for polygon in self.polygons:
            path = Path(polygon.vertices)
            mask &= ~path.contains_points(points)
        
        # Вычисляем поле только для точек вне полигонов
        Ex = np.zeros(X.shape).ravel()
        Ey = np.zeros(Y.shape).ravel()
        valid_points = points[mask]
        if len(valid_points) > 0:
            field_x, field_y = zip(*[self.calculate_total_field(p[0], p[1]) for p in valid_points])
            Ex[mask] = field_x
            Ey[mask] = field_y
        
        Ex = Ex.reshape(X.shape)
        Ey = Ey.reshape(Y.shape)
        
        self.ax.clear()
        E_magnitude = np.sqrt(Ex**2 + Ey**2)
        
        # Оптимизируем параметры streamplot
        self.ax.streamplot(X, Y, Ex, Ey,
                          color=E_magnitude,
                          cmap='viridis',
                          density=1.2,  # уменьшена плотность
                          linewidth=1.2,
                          arrowsize=0.8,
                          arrowstyle='->',
                          broken_streamlines=False,
                          minlength=0.1)

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
        dt = 0.001
        
        # Рисуем поле один раз
        self.visualize_field()
        
        # Удаляем начальную точку пробного заряда
        if len(self.ax.lines) > 0:
            self.ax.lines[-1].remove()
        
        while self.animation_running:
            test_body_mass = self.test_body.mass/100;
            Ex, Ey = self.calculate_total_field(self.test_body.x, self.test_body.y)
            scale_factor = 1e-8
            Fx = Ex * self.test_body.charge * scale_factor
            Fy = Ey * self.test_body.charge * scale_factor

            self.test_body.update_position(Fx, Fy, dt)

            # Обновляем только положение тестового заряда
            self.ax.plot(self.test_body.x, self.test_body.y, 'o', color='green', markersize=8)
            self.canvas.draw()
            
            # Удаляем старую точку
            self.ax.lines[-1].remove()
            # Проверяем, не вышел ли заряд за пределы экрана
            if (self.test_body.x < 0 or self.test_body.x > 1 or 
                self.test_body.y < 0 or self.test_body.y > 1):
                self.animation_running = False
                return
            self.root.update()
            # time.sleep(0.0005)  # Замедляем анимацию для лучшей видимости


if __name__ == "__main__":
    root = Tk()
    app = ElectricFieldApp(root)
    root.mainloop()
