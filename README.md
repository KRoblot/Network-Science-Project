

# **Network Science Project – Epidemic Spreading on Signed Networks**  

## **📌 Overview**  
This project aims to **critically verify** the claims presented in the paper:  
📄 *"Epidemic Spreading on Evolving Signed Networks"* (Phys. Rev. E 95, 022314, 2017).  

The paper explores how **social relationships evolve in response to an epidemic**, using a **signed network model** where links can be either **friendly (+1)** or **hostile (-1)**. Our goal is to **reproduce the key results, analyze their validity, and test the robustness of their conclusions** under different conditions with a critical point of view.  

## **🎯 Objectives** 
✅ Reproduce the **figures and key numerical results** from the paper.  
✅ Test the **impact of different network structures** (beyond complete graphs).  
✅ Analyze the **role of key parameters** (infection rate \( \alpha \), initial infected proportion \( \rho_0 \), etc.).  
✅ Investigate **jammed states** and the **energy function's validity**.  
✅ Identify **potential limitations** and propose improvements to the model.  

## **📂 Repository Structure**  
```
📦 Network-Science-Project
├── 📄requirements.txt
├── 📜 README.md                    # Project documentation  
├── 📂 __pycache__/                  # Compiled Python cache files  
├── 📜 main.py                       # Main script to run simulations  
├── 📜 calculate.py                  # Core functions for network calculations  
├── 📜 verify_jammed_state.py        # Script to verify the existence of jammed states  
├── 📜 density_of_jammed_state.py    # Computes the density of jammed states  
├── 📜 density_of_jammed_state_alt.py# Alternative method for jammed state density  
├── 📜 evolution_density_infected_nodes.py  # Tracks infected nodes over time  
├── 📜 fig8.py                        # Script to reproduce Figure 8 from the paper  
├── 📜 fig10.py                       # Script to reproduce Figure 10 from the paper  
├── 📜 'Figure 8.py'                  # Old script for Figure 8  
├── 📂 figures/                       # Directory for generated figures  
│   ├── density_of_jammed_state.png   # Generated density of jammed states  
│   ├── evolution_of_infected_nodes.png # Evolution of infected nodes  
│   ├── 3dfigure.png                  # 3D visualization  
│   ├── graph.png                      # Graph representation  

## **📊 Methodology**  
1. **Network Initialization**  
   - Create a signed **Erdős–Rényi** or **complete graph** where edges are assigned friendly (+1) or hostile (-1) signs.  
   - Define **initial infected nodes** based on proportion \( \rho_0 \).  

2. **Epidemic Dynamics (SI Model on Signed Networks)**  
   - The disease spreads **only through friendly edges** with probability \( \alpha \).  
   - The network evolves dynamically:  
     - Friendly edges **may become hostile** if infection occurs.  
     - Nodes with the same state **reinforce their friendships**.  

3. **Monte Carlo Simulations & Energy Function**  
   - Compute **structural balance** and **jammed states** using an energy function.  
   - Investigate **alternative definitions** for stopping conditions (e.g., network fragmentation).  

4. **Critical Analysis**  
   - Compare **simulation results with the paper’s predictions**.  
   - Identify **discrepancies, biases, and missing real-world factors**.  
   - Test **alternative epidemic models** (e.g., adding recovery in SIR).  

## **📌 How to Run the Project**  
### **1️⃣ Install Dependencies**  
```bash
pip install -r requirements.txt
```
### **2️⃣ Run the Main Simulations**  
```bash
python src/density_of_jammed_state.py
python src/density_of_jammed_state_alt.py
python src/evolution_density_infected_nodes.py
python src/fig8.py
python src/fig10.py
```


## **📈 Expected Outcomes & Open Questions**  
🔹 Does **social tension minimization** actually lead to jammed states as described?  
🔹 How does **network structure** (e.g., scale-free vs. complete) affect epidemic outcomes?  
🔹 Can **small rule modifications** (e.g., recovery, delayed reactions) significantly alter conclusions?  
🔹 Are **Monte Carlo update rules (p = 0.5 for neutral moves)** justified?  
🔹 Can we **bridge the gap** between this theoretical model and real-world epidemic data?  

## **🛠 Contributors**  
- **Koské Roblot and Debbah Adam** - Main development & analysis  
- **Sirot Marine and Wintrebert Baptiste** - Contributions in simulations, theory, and review  
- **Bongiorno Christian** - Supervisor

This project is licensed under the **MIT License**.  

MIT License  

Copyright (c) 2025 Koské Roblot, Debbah Adam, Sirot Marine, Wintrebert Baptiste 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.  
