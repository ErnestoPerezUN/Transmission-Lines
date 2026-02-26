import tkinter as tk
from tkinter import filedialog, messagebox
import os

# Example scripts
def script_one(file_path):
    messagebox.showinfo("Script 1", f"Running Script 1 on:\n{file_path}")

def script_two(file_path):
    messagebox.showinfo("Script 2", f"Running Script 2 on:\n{file_path}")

class FileApp:
    def __init__(self, root):
        self.root = root
        self.root.title("File Selector & Script Runner")
        self.root.geometry("400x200")

        self.file_path = None

        # Select File button
        self.select_button = tk.Button(root, text="Select File", command=self.select_file)
        self.select_button.pack(pady=10)

        # Label to show chosen file
        self.label = tk.Label(root, text="No file selected")
        self.label.pack(pady=5)

        # Script buttons
        self.btn_script1 = tk.Button(root, text="Run Script 1", command=self.run_script1)
        self.btn_script1.pack(pady=5)

        self.btn_script2 = tk.Button(root, text="Run Script 2", command=self.run_script2)
        self.btn_script2.pack(pady=5)

    def select_file(self):
        file_path = filedialog.askopenfilename()
        if file_path:
            self.file_path = file_path
            self.label.config(text=os.path.basename(file_path))

    def run_script1(self):
        if self.file_path:
            script_one(self.file_path)
        else:
            messagebox.showwarning("Warning", "Please select a file first.")

    def run_script2(self):
        if self.file_path:
            script_two(self.file_path)
        else:
            messagebox.showwarning("Warning", "Please select a file first.")

if __name__ == "__main__":
    root = tk.Tk()
    app = FileApp(root)
    root.mainloop()