import pandas as pd
import matplotlib.pyplot as plt
import pathlib

rel_path = pathlib.Path(__file__).parents[0] # relative_path

'''параметры настройки'''
CSV_FILE = rel_path / '..\\build-vscode\\tests\\Debug\data.csv'           # Путь к CSV файлу
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
plt.plot(list(df['r']), list(df['Tref']), linestyle=LINE_STYLE, color='blue', label='reference')
plt.plot(list(df['r']), list(df['Tcalc']), linestyle=LINE_STYLE, color='red', marker=MARKER, label='calculated')

plt.title(PLOT_TITLE)
plt.xlabel(X_LABEL)
plt.ylabel(Y_LABEL)

if SHOW_LEGEND:
    plt.legend()

plt.grid(True)
plt.tight_layout()
plt.show()
