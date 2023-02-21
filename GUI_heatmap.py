import tkinter as tk
from tkinter import ttk
from tkinter import *
from tkinter import filedialog
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit import Chem

class Application(tk.Frame):
    def __init__(self,master):
        super().__init__(master)
        self.pack()
        self.master.geometry("1000x800")
        self.master.title("Python GUI Tool")
        self.create_widgets()

    # ウィジェットの作成
    def create_widgets(self):
        # Frame1の作成
        self.frame1 = ttk.Frame(self, padding=10)
        self.frame1.grid()

        # 「ファイル」ラベルの作成
        self.s = StringVar()
        self.s.set('ファイル名：')
        self.label1 = ttk.Label(self.frame1, textvariable=self.s)
        self.label1.grid(row=0, column=0)

        # 参照ファイルのパスを表示するテキストボックスの作成
        self.file_path = StringVar()
        self.filepath_entry = ttk.Entry(self.frame1, textvariable=self.file_path, width=50)
        self.filepath_entry.grid(row=0, column=1)

        # 参照ボタンの作成
        self.refer_button = ttk.Button(self.frame1, text=u'参照', command=self.click_refer_button)
        self.refer_button.grid(row=0, column=2)

        # Frame2の作成
        self.frame2 = ttk.Frame(self, padding=10)
        self.frame2.grid()

        # テキスト出力ボタンの作成
        self.export_button = ttk.Button(self.frame2, text='ファイルの中身を出力', command=self.click_export_button, width=15)
        self.export_button.grid(row=0, column=0)
  
        # テキスト出力ボックスの作成
        self.textboxname = StringVar()
        self.textboxname.set('\n\n出力内容 ')
        self.label3 = ttk.Label(self.frame2, textvariable=self.textboxname)
        self.label3.grid(row=1, column=0)
        self.textBox = Text(self.frame2, width=80)
        self.textBox.grid(row=2, column=0)

        #Frame3の作成
        self.frame3 = ttk.Frame(self, padding=10)
        self.frame3.grid()

        #runボタンの作成
        self.run_button = ttk.Button(self.frame3, text='実行', command=lambda:[self.delete_frames(), self.graph()], width=10)
        self.run_button.grid(row=9, column=0)

    # ファイルの参照処理
    def click_refer_button(self):
        self.fTyp = [("","*")]
        self.iDir = os.path.abspath(os.path.dirname(__file__))
        self.filepath = filedialog.askopenfilename(filetypes = self.fTyp, initialdir = self.iDir)
        self.file_path.set(self.filepath)

    # 出力処理
    def click_export_button(self):
        self.path = self.file_path.get()
        self.f = open(self.path, encoding="utf-8")
        self.text_data = self.f.read()
        self.textBox.insert(END, self.text_data)

    def graph(self):
        self.frame4 = ttk.Frame(self)
        self.frame4.grid()
        suppl = Chem.SDMolSupplier(self.filepath)
        mols = [mol for mol in suppl if mol is not None]
        fps = [AllChem.GetMorganFingerprint(mol, 2, useFeatures=True) for mol in mols]
        sim_matrix = [DataStructs.BulkTanimotoSimilarity(fp, fps) for fp in fps]

        fig, ax = plt.subplots()
        im = ax.imshow(sim_matrix, cmap='bwr')
        plt.colorbar(im)

        #Figureを埋め込み
        canvas = FigureCanvasTkAgg(fig, self.frame4)
        canvas.get_tk_widget().pack()
        #ツールバーを表示
        toolbar=NavigationToolbar2Tk(canvas, self.frame4)

    def delete_frames(self):
        self.frame1.destroy()
        self.frame2.destroy()
        self.frame3.destroy()

if __name__ == "__main__":
    root = tk.Tk()
    app = Application(master = root)
    app.mainloop()