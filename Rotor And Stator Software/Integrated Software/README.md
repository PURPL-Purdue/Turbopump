# Setup Instructions

```bash
cd "./Turbopump/Rotor And Stator Software/Integrated Software"
conda create --name turbopump python=3.12
conda activate turbopump
pip install -r requirements.txt
cd turbo-design-main
pip install .
cd ..
python main.py
```