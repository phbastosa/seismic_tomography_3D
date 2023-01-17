import os
import tkinter
import customtkinter

from PIL import Image

customtkinter.set_appearance_mode("System")  # Modes: "System" (standard), "Dark", "Light"
customtkinter.set_default_color_theme("dark-blue")  # Themes: "blue" (standard), "green", "dark-blue"

class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        # configure window
        self.width = 1500
        self.height = 840

        self.title("3D OBN First Arrival Tomography")
        self.geometry(f"{self.width}x{self.height}+200+100")  
        self.resizable(False, False)

        self.menu_frame = customtkinter.CTkFrame(self, width = 200, corner_radius=0)
        self.menu_frame.place(x = 0, y = 0, relheight = 1)

        image_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "img/")
        self.logo_image = customtkinter.CTkImage(Image.open(os.path.join(image_path, "logo.png")), size=(180, 130))

        self.logo_frame = customtkinter.CTkFrame(self.menu_frame, height = 130, width = 180, corner_radius=0) 
        self.logo_frame.place(x = 10, y = 10)

        self.img_frame_label = customtkinter.CTkLabel(self.logo_frame, text = " ", image=self.logo_image)
        self.img_frame_label.place(x = 0, y = 0)

        self.name_frame = customtkinter.CTkFrame(self.logo_frame, fg_color = ("#EEEEEE","#191A19"), height = 30, width = 180, corner_radius=0)    
        self.name_frame.place(x = 0, y = 100)

        self.name_label = customtkinter.CTkLabel(self.name_frame, text="3D HeadTomoGrid", font=customtkinter.CTkFont(size=15))
        self.name_label.place(relx = 0.5, rely = 0.5, anchor = "center")

        # Application buttons
        self.model_buttom = customtkinter.CTkButton(self.menu_frame, height = 40, width = 160, text = "Model", font=customtkinter.CTkFont(size=20, weight="normal")) #, command=self.model_button_event)
        self.model_buttom.place(x = 20, rely = 200/self.height)

        self.geometry_buttom = customtkinter.CTkButton(self.menu_frame, height = 40, width = 160, text = "Geometry", font=customtkinter.CTkFont(size=20, weight="normal"))#, command=self.geometry_button_event)
        self.geometry_buttom.place(x = 20, rely = 260/self.height)

        self.eikonal_buttom = customtkinter.CTkButton(self.menu_frame, height = 40, width = 160, text = "Eikonal", font=customtkinter.CTkFont(size=20, weight="normal"))#, command=self.source_button_event)
        self.eikonal_buttom.place(x = 20, rely = 320/self.height)

        self.tomography_buttom = customtkinter.CTkButton(self.menu_frame, height = 40, width = 160, text = "Tomography", font=customtkinter.CTkFont(size=20, weight="normal"))#, command=self.tomography_button_event)
        self.tomography_buttom.place(x = 20, rely = 380/self.height)

        self.appearance_mode_label = customtkinter.CTkLabel(self.menu_frame, text="Appearance Mode", font=customtkinter.CTkFont(size=17))
        self.appearance_mode_label.place(relx = 0.5, rely = 0.9, anchor = "center")
        self.appearance_mode_optionemenu = customtkinter.CTkOptionMenu(self.menu_frame, height = 30, width = 160, font=customtkinter.CTkFont(size=15), values=["Light", "Dark", "System"],command=self.change_appearance_mode_event)
        self.appearance_mode_optionemenu.place(relx = 0.5, rely = 0.95, anchor = "center")
        self.appearance_mode_optionemenu.set("System")
    
    def change_appearance_mode_event(self, new_appearance_mode: str):
        customtkinter.set_appearance_mode(new_appearance_mode)

#     def main_frame_clear(self):
#         for widget in self.main_frame.winfo_children():
#             widget.destroy()

#         self.main_frame.place_forget()
#         self.main_frame.place(relx = 220/self.width + 0.02, rely = 0.015, relwidth = 0.787, relheight = 0.97)    


# #-----------------------------------------------------------------------------------------
# # Velocity model -------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------

#     def model_button_event(self):
#         print("model_button click")
#         self.main_frame_clear()

#         self.model_dimension_frame = customtkinter.CTkFrame(self.main_frame) 
#         self.model_dimension_frame.place(relx = 0.02, rely = 0.035, relwidth = 0.47, relheight = 0.45)    

#         self.model_layerCake_frame = customtkinter.CTkFrame(self.main_frame) 
#         self.model_layerCake_frame.place(relx = 0.51, rely = 0.035, relwidth = 0.47, relheight = 0.45)    

#         self.model_gradient_frame = customtkinter.CTkFrame(self.main_frame) 
#         self.model_gradient_frame.place(relx = 0.02, rely = 0.515, relwidth = 0.47, relheight = 0.45)    

#         self.model_visualize_frame = customtkinter.CTkFrame(self.main_frame) 
#         self.model_visualize_frame.place(relx = 0.51, rely = 0.515, relwidth = 0.47, relheight = 0.45)    

# #-----------------------------------------------------------------------------------------
# # Geometry -------------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------

#     def geometry_button_event(self):
#         print("geometry_button click")
#         self.main_frame_clear()


# #-----------------------------------------------------------------------------------------
# # Source ---------------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------

#     def source_button_event(self):
#         print("source_button click")
#         self.main_frame_clear()

# #-----------------------------------------------------------------------------------------
# # Acoustic -------------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------

#     def acoustic_button_event(self):
#         print("acoustic_button click")
#         self.main_frame_clear()

# #-----------------------------------------------------------------------------------------
# # Auto-picking ---------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------

#     def autoPick_button_event(self):
#         print("autoPick_button click")
#         self.main_frame_clear()

# #-----------------------------------------------------------------------------------------
# # Eikonal --------------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------

#     def eikonal_button_event(self):
#         print("eikonal_button click")
#         self.main_frame_clear()

# #-----------------------------------------------------------------------------------------
# # Tomography -----------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------

#     def tomography_button_event(self):
#         print("tomography_button click")
#         self.main_frame_clear()


#-----------------------------------------------------------------------------------------
# Main function --------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

if __name__ == "__main__":
    app = App()
    app.mainloop()
