from tkinter import Tk
from src.gui.electric_field_app import ElectricFieldApp

if __name__ == "__main__":
    root = Tk()
    app = ElectricFieldApp(root)
    root.mainloop()