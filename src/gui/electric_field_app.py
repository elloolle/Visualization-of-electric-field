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
    def __init__(self, root):
        self.root = root
        self.root.title("Визуализация электрического поля")

        self.charges = []
        self.polygons = []
        self.current_polygon = None
        self.field_density = 20
        self.test_body = None
        self.placing_test_charge = False
        self.animation_running = False

        self._create_control_elements()
        self._setup_matplotlib()
        self._add_initial_charge()

    def _create_control_elements(self):
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
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        self.canvas.mpl_connect("button_press_event", self.add_point)
        self.canvas.mpl_connect("key_press_event", self.close_polygon)

        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)

    def _add_initial_charge(self):
        initial_charge = Charge(0.5, 0.5, 1)
        self.charges.append(initial_charge)
        self.ax.plot(0.5, 0.5, 'o', color='red')
        self.canvas.draw()

    def start_placing_test_charge(self):
        self.placing_test_charge = True

    def add_point(self, event):
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
        charge_value = simpledialog.askfloat("Величина заряда", "Введите величину заряда:", parent=self.root)
        if charge_value is not None:
            charge = Charge(x, y, charge_value)
            self.charges.append(charge)
            color = 'red' if charge_value > 0 else 'blue'
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
        self.field_density = int(self.density_scale.get() * 20)
        
        x_coords = np.linspace(0, 1, self.field_density)
        y_coords = np.linspace(0, 1, self.field_density)
        X, Y = np.meshgrid(x_coords, y_coords)
        points = np.column_stack((X.ravel(), Y.ravel()))
        
        mask = np.ones(len(points), dtype=bool)
        for polygon in self.polygons:
            path = Path(polygon.vertices)
            mask &= ~path.contains_points(points)
        
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
        self.ax.clear()
        field_magnitude = np.sqrt(Ex**2 + Ey**2)
        
        self.ax.streamplot(X, Y, Ex, Ey,
                          color=field_magnitude,
                          cmap='viridis',
                          density=1,
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
        TIME_STEP = 0.001
        FORCE_SCALE = 1e-8
        
        self.visualize_field()
        
        if len(self.ax.lines) > 0:
            self.ax.lines[-1].remove()
        
        while self.animation_running:
            scaled_mass = self.test_body.mass/100
            field_x, field_y = self.calculate_total_field(self.test_body.x, self.test_body.y)
            force_x = field_x * self.test_body.charge * FORCE_SCALE
            force_y = field_y * self.test_body.charge * FORCE_SCALE

            self.test_body.update_position(force_x, force_y, TIME_STEP)

            self.ax.plot(self.test_body.x, self.test_body.y, 'o', color='green', markersize=8)
            self.canvas.draw()
            
            self.ax.lines[-1].remove()
            
            if not (0 <= self.test_body.x <= 1 and 0 <= self.test_body.y <= 1):
                self.animation_running = False
                return
                
            self.root.update() 