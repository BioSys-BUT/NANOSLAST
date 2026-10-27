import customtkinter as ctk
from .search_part import create_search_part
from .plot_part import create_plot_part
from .export_part import create_export_part
from .plot_specific_part import create_plot_specific_signal_part

ctk.set_appearance_mode("System")
ctk.set_default_color_theme("blue")

class App(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("NANOSLAST")
        self.geometry("1000x500")

        self.tabview = ctk.CTkTabview(self,
                                      width=980,
                                      height=640,
                                      corner_radius=10)
        self.tabview.pack(expand=True, fill="both", padx=20, pady=20)

        self.tabview.add("Search")
        self.tabview.add("Export")
        self.tabview.add("Plot")
        # self.tabview.add("Plot Specific Signal")

        create_search_part(self.tabview.tab("Search"))
        create_export_part(self.tabview.tab("Export"))
        create_plot_part(self.tabview.tab("Plot"))
        # create_plot_specific_signal_part(self.tabview.tab("Plot Specific Signal"))  # Přidání nové části


        self.tabview._segmented_button.configure(
            font=ctk.CTkFont(size=16, weight="bold")
        )

        for btn in self.tabview._segmented_button._buttons_dict.values():
            btn.configure(width=120)
