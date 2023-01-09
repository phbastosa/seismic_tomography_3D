import tkinter
import tkinter.messagebox
import customtkinter

customtkinter.set_appearance_mode("System")  # Modes: "System" (standard), "Dark", "Light"
customtkinter.set_default_color_theme("blue")  # Themes: "blue" (standard), "green", "dark-blue"

class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        # configure window
        self.width = 1200
        self.height = 600
        self.title("3D seismic tomography")
        self.geometry(f"{self.width}x{self.height}")  
        self.resizable(False, False)
            
        self.title_frame = customtkinter.CTkFrame(self, fg_color = "transparent", height = 50, width = 220)
        self.title_frame.place(relx = 0.01, rely = 0.01, relwidth = 220/self.width, relheight = 50/self.height)
        
        self.applications_label = customtkinter.CTkLabel(self.title_frame, text="Applications", font=customtkinter.CTkFont(size=30))
        self.applications_label.place(relx = 0.5, rely = 0.5, anchor = "center")

        self.line1_frame = customtkinter.CTkFrame(self, height = 5, width = 220)
        self.line1_frame.place(relx = 0.01, rely = 0.1, relwidth = 220/self.width, relheight = 5/self.height)

        # Application buttons
        self.applications_frame = customtkinter.CTkFrame(self, height = 350, width = 220)
        self.applications_frame.place(relx = 0.01, rely = 0.12, relwidth = 220/self.width, relheight = 400/self.height)

        self.model_buttom = customtkinter.CTkButton(self.applications_frame, height = 50, width = 200, text = "Velocity model", font=customtkinter.CTkFont(size=20), command=self.model_button_event)
        self.model_buttom.place(relx = 0.5, rely = 0.5 * 40/self.applications_frame.winfo_reqheight() + 0.02, relwidth = 200/self.applications_frame.winfo_reqwidth(), relheight = 40/self.applications_frame.winfo_reqheight(), anchor = "center")

        self.geometry_buttom = customtkinter.CTkButton(self.applications_frame, height = 50, width = 200, text = "OBN Geometry", font=customtkinter.CTkFont(size=20), command=self.geometry_button_event)
        self.geometry_buttom.place(relx = 0.5, rely = 0.5 * 40/self.applications_frame.winfo_reqheight() + 0.16, relwidth = 200/self.applications_frame.winfo_reqwidth(), relheight = 40/self.applications_frame.winfo_reqheight(), anchor = "center")

        self.source_buttom = customtkinter.CTkButton(self.applications_frame, height = 50, width = 200, text = "Source signature", font=customtkinter.CTkFont(size=20), command=self.source_button_event)
        self.source_buttom.place(relx = 0.5, rely = 0.5 * 40/self.applications_frame.winfo_reqheight() + 0.30, relwidth = 200/self.applications_frame.winfo_reqwidth(), relheight = 40/self.applications_frame.winfo_reqheight(), anchor = "center")

        self.acoustic_buttom = customtkinter.CTkButton(self.applications_frame, height = 50, width = 200, text = "Acoustic wave", font=customtkinter.CTkFont(size=20), command=self.acoustic_button_event)
        self.acoustic_buttom.place(relx = 0.5, rely = 0.5 * 40/self.applications_frame.winfo_reqheight() + 0.44, relwidth = 200/self.applications_frame.winfo_reqwidth(), relheight = 40/self.applications_frame.winfo_reqheight(), anchor = "center")

        self.autoPick_buttom = customtkinter.CTkButton(self.applications_frame, height = 50, width = 200, text = "Automatic picking", font=customtkinter.CTkFont(size=20), command=self.autoPick_button_event)
        self.autoPick_buttom.place(relx = 0.5, rely = 0.5 * 40/self.applications_frame.winfo_reqheight() + 0.58, relwidth = 200/self.applications_frame.winfo_reqwidth(), relheight = 40/self.applications_frame.winfo_reqheight(), anchor = "center")

        self.eikonal_buttom = customtkinter.CTkButton(self.applications_frame, height = 50, width = 200, text = "Eikonal solvers", font=customtkinter.CTkFont(size=20), command=self.eikonal_button_event)
        self.eikonal_buttom.place(relx = 0.5, rely = 0.5 * 40/self.applications_frame.winfo_reqheight() + 0.72, relwidth = 200/self.applications_frame.winfo_reqwidth(), relheight = 40/self.applications_frame.winfo_reqheight(), anchor = "center")

        self.tomography_buttom = customtkinter.CTkButton(self.applications_frame, height = 50, width = 200, text = "Seismic tomography", font=customtkinter.CTkFont(size=20), command=self.tomography_button_event)
        self.tomography_buttom.place(relx = 0.5, rely = 0.5 * 40/self.applications_frame.winfo_reqheight() + 0.86, relwidth = 200/self.applications_frame.winfo_reqwidth(), relheight = 40/self.applications_frame.winfo_reqheight(), anchor = "center")

        self.line2_frame = customtkinter.CTkFrame(self, height = 5, width = 220)
        self.line2_frame.place(relx = 0.01, rely = 0.8, relwidth = 220/self.width, relheight = 5/self.height)

        self.config_frame = customtkinter.CTkFrame(self)
        self.config_frame.place(relx = 0.01, rely = 0.82, relwidth = 220/self.width, relheight = 100/self.height)

        self.appearance_mode_label = customtkinter.CTkLabel(self.config_frame, text="Appearance Mode", font=customtkinter.CTkFont(size=20))
        self.appearance_mode_label.place(relx = 0.5, rely = 0.2, anchor = "center")
        self.appearance_mode_optionemenu = customtkinter.CTkOptionMenu(self.config_frame, height = 40, width = 200, font=customtkinter.CTkFont(size=20), values=["Light", "Dark", "System"],command=self.change_appearance_mode_event)
        self.appearance_mode_optionemenu.place(relx = 0.5, rely = 0.65, anchor = "center")
        self.appearance_mode_optionemenu.set("Dark")

        self.main_frame = customtkinter.CTkFrame(self)

    def main_frame_clear(self):
        for widget in self.main_frame.winfo_children():
            widget.destroy()

        self.main_frame.place_forget()
        self.main_frame.place(relx = 220/self.width + 0.02, rely = 0.015, relwidth = 0.787, relheight = 0.97)    

    def change_appearance_mode_event(self, new_appearance_mode: str):
        customtkinter.set_appearance_mode(new_appearance_mode)

#-----------------------------------------------------------------------------------------
# Velocity model -------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

    def model_button_event(self):
        print("model_button click")
        self.main_frame_clear()

        self.model_dimension_frame = customtkinter.CTkFrame(self.main_frame) 
        self.model_dimension_frame.place(relx = 0.02, rely = 0.035, relwidth = 0.47, relheight = 0.45)    

        self.model_layerCake_frame = customtkinter.CTkFrame(self.main_frame) 
        self.model_layerCake_frame.place(relx = 0.51, rely = 0.035, relwidth = 0.47, relheight = 0.45)    

        self.model_gradient_frame = customtkinter.CTkFrame(self.main_frame) 
        self.model_gradient_frame.place(relx = 0.02, rely = 0.515, relwidth = 0.47, relheight = 0.45)    

        self.model_visualize_frame = customtkinter.CTkFrame(self.main_frame) 
        self.model_visualize_frame.place(relx = 0.51, rely = 0.515, relwidth = 0.47, relheight = 0.45)    

#-----------------------------------------------------------------------------------------
# Geometry -------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

    def geometry_button_event(self):
        print("geometry_button click")
        self.main_frame_clear()


#-----------------------------------------------------------------------------------------
# Source ---------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

    def source_button_event(self):
        print("source_button click")
        self.main_frame_clear()

#-----------------------------------------------------------------------------------------
# Acoustic -------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

    def acoustic_button_event(self):
        print("acoustic_button click")
        self.main_frame_clear()

#-----------------------------------------------------------------------------------------
# Auto-picking ---------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

    def autoPick_button_event(self):
        print("autoPick_button click")
        self.main_frame_clear()

#-----------------------------------------------------------------------------------------
# Eikonal --------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

    def eikonal_button_event(self):
        print("eikonal_button click")
        self.main_frame_clear()

#-----------------------------------------------------------------------------------------
# Tomography -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

    def tomography_button_event(self):
        print("tomography_button click")
        self.main_frame_clear()


#-----------------------------------------------------------------------------------------
# Main function --------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

if __name__ == "__main__":
    app = App()
    app.mainloop()