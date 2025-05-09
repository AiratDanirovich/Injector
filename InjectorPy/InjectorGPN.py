import os
import json

import InjectorGPN

import pandas as pd

# help(InjectorGPN)  # посмотреть, что доступно

class Solver:
    def __init__(self, config=None):
        if not config:
            self.load_config()
        else:
            self.config = config
        t1 = self.config["t_start"] + 2e6
        time_step = (t1 - self.config["t_start"]) / self.config["time_steps_nmbr"]
        self.config["time_intervals"] = self.config["time_steps_nmbr"]*[time_step]
        self.result = {}

    def load_config(self):
        with open(os.path.join('input', 'config.json'), 'r') as fin:
            self.config = json.load(fin)

    def main_calc(self):

        wrapper = InjectorGPN.Wrapper(
            self.config["q"],
            self.config["density"],
            self.config["capacity"],
            self.config["viscosity"],
            self.config["rMin"],
            self.config["rMax"],
            self.config["r_nodes_nmbr"],
            self.config["zTop"],
            self.config["thickness"],
            self.config["conductivity"],
            self.config["porosity"],
            self.config["is_permeable"],
            self.config["solid_density"],
            self.config["solid_specific_heatcapacity"],
            self.config["t_start"],
            self.config["time_intervals"]
        )

        # Получаем и выводим результаты
        times = wrapper.get_times()
        temps = wrapper.get_temps()
        self.result = {"temps_{}".format(i): t_row for i, t_row in enumerate(temps)}

        # print("Временные шаги:", times)
        # print("Температурные профили (по времени и радиусу):")

        self.export_result()

    def export_result(self):
        pd.DataFrame(self.result).to_excel(os.path.join('output', 'result.xlsx'), index=False)

def main ():
    sc = Solver()
    sc.main_calc()

if __name__ == '__main__':
    main()