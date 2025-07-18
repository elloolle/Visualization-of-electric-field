# Визуализация линий напряженности электрического поля системы зарядов

## Описание проекта
Эта программа позволяет пользователю отмечать заряды и замкнутые области (многоугольники) с равномерно распределенными зарядами, задавать их величину и визуализировать линии напряженности электрического поля системы зарядов. Пользователь может ввести систему зарядов, а затем, с помощью кнопки "Старт визуализации", программа отобразит линии напряженности электрического поля.

## Реализуемый функционал
1. **Запуск программы**:
   - Открывается окно с кнопкой "Старт визуализации" и ползунком "Частота линий" (от 0 до 1), который регулирует количество линий напряженности на экране.
   
2. **Добавление заряда**:
   - При клике мышью появляется окно для ввода величины заряда. В окне есть поле ввода и две кнопки: "Применить" и "Отмена".
   - После ввода значения заряда и нажатия "Применить" на экране появляется точка синего или красного цвета (в зависимости от знака заряда) с отображением значения заряда.

3. **Замкнутые области**:
   - Пользователь может отметить замкнутую область (многоугольник), которая заполняется синим или красным цветом, обозначая равномерно распределенный по площади заряд внутри этой области.
   
4. **Визуализация линий напряженности**:
   - После добавления зарядов и многоугольников пользователь нажимает кнопку "Старт визуализации", и программа отображает линии напряженности электрического поля.
   
5. **Движение заряженных частиц**:
   - Пользователь может выбрать точку, задать ей заряд и массу. При нажатии кнопки "Старт" эта частица начинает двигаться под действием сил электрического поля.

## Архитектура
Проект включает пять основных классов: `Charge`, `Vector`, `Point`, `Polygon` и `Body`.

### Класс `Charge`
- **Атрибуты**:
  - `x`, `y` — координаты заряда.
  - `value` — величина заряда.

- **Методы**:
  - `electric_field(x1, y1)` — возвращает вектор напряженности электрического поля заряда в точке `(x1, y1)`.

### Класс `Vector`
- **Атрибуты**:
  - `x1`, `y1`, `x2`, `y2` — координаты начала и конца вектора.

- **Методы**:
  - `sum(v1, v2, e)` — возвращает вектор, равный сумме `v1 + v2`, предполагая, что конец `v1` совпадает с началом `v2` с погрешностью `e`.
  - `end()` — возвращает точку, соответствующую концу вектора.

### Класс `Point`
- **Атрибуты**:
  - `x`, `y` — координаты точки.

### Класс `Polygon`
Этот класс представляет замкнутую область, где заряд распределен равномерно.

- **Атрибуты**:
  - `vertices` — список объектов `Point`, представляющих вершины многоугольника.
  - `charge_density` — плотность заряда на единицу площади внутри многоугольника.

- **Методы**:
  - `add_vertex(x, y)` — добавляет новую вершину с координатами `(x, y)` в список `vertices`.
  - `is_closed()` — возвращает `True`, если первая и последняя точки совпадают (область замкнута).
  - `fill_color()` — возвращает синий или красный цвет для заливки в зависимости от знака `charge_density`.
  - `electric_field(x1, y1)` — вычисляет напряженность электрического поля многоугольника в точке `(x1, y1)`, предполагая равномерное распределение заряда.
  - `break_into_charges(width, length)` — разбивает многоугольник на множество точечных зарядов с заданной шириной и длиной. Возвращает массив объектов типа `Charge`, которые могут быть рассмотрены как точечные заряды.

### Класс `Body`
Этот класс представляет физическое тело, которое будет иметь координаты, скорость и ускорение.

- **Атрибуты**:
  - `position` — вектор, представляющий текущие координаты тела (например, `(x, y)`).
  - `velocity` — вектор, представляющий скорость тела, указывающий направление и величину движения.
  - `acceleration` — вектор, представляющий ускорение тела, указывающий изменения скорости со временем.

- **Методы**:
  - `update(dt)` — обновляет положение и скорость тела на основе времени `dt`. Этот метод будет использоваться для расчета новых координат тела в зависимости от его текущей скорости и ускорения.
  - `apply_force(force_vector)` — применяет силу (вектор) к телу, что изменяет его ускорение. Этот метод позволит взаимодействовать с телом, например, под воздействием электрических сил.
