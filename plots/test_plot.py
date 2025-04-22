import pandas as pd
import matplotlib.pyplot as plt

'''параметры настройки'''
CSV_FILE = 'D:\Lessons\Grants\\2025\GPN\Injector\plots\data.csv'           # Путь к CSV файлу
DELIMETER = ';'                 # Разделитель
PLOT_TITLE = 'График'           # Название графика
X_LABEL = 'X'                   # Подпись оси X
Y_LABEL = 'Y'                   # Подпись оси Y
SHOW_LEGEND = True              # Показывать легенду
LEGEND_LABEL = 'Значения Y'     # Подпись в легенде
LINE_STYLE = '-'                # Стиль линии (например '-', '--', '-.', ':')
MARKER = 'o'                    # Маркер на точках (например 'o', 's', '^', '')
COLOR = 'blue'                  # Цвет линии

# импорт csv
df = pd.read_csv(CSV_FILE, delimiter=DELIMETER, header=0)

# нарисовать график
plt.figure(figsize=(10, 6))
plt.plot(list(df['x']), list(df['y']), linestyle=LINE_STYLE, marker=MARKER, color=COLOR, label=LEGEND_LABEL)

plt.title(PLOT_TITLE)
plt.xlabel(X_LABEL)
plt.ylabel(Y_LABEL)

if SHOW_LEGEND:
    plt.legend()

plt.grid(True)
plt.tight_layout()
plt.show()
