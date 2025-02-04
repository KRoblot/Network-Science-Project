

# **Network Science Project – Epidemic Spreading on Signed Networks**  

## **📌 Overview**  
This project aims to **critically verify** the claims presented in the paper:  
📄 *"Epidemic Spreading on Evolving Signed Networks"* (Phys. Rev. E 95, 022314, 2017).  

The paper explores how **social relationships evolve in response to an epidemic**, using a **signed network model** where links can be either **friendly (+1)** or **hostile (-1)**. Our goal is to **reproduce the key results, analyze their validity, and test the robustness of their conclusions** under different conditions.  

## **🎯 Objectives**  
✅ Reproduce the **figures and key numerical results** from the paper.  
✅ Test the **impact of different network structures** (beyond complete graphs).  
✅ Analyze the **role of key parameters** (infection rate \( \alpha \), initial infected proportion \( \rho_0 \), etc.).  
✅ Investigate **jammed states** and the **energy function's validity**.  
✅ Identify **potential limitations** and propose improvements to the model.  

## **📂 Repository Structure**  
```
📦 Network-Science-Project  
├── 📜 README.md                # Project documentation  
├── 📂 src                      # Main source code  
│   ├── epidemic_model.py       # Implementation of the SI model on signed networks  
│   ├── monte_carlo_sim.py      # Monte Carlo simulations for network evolution  
│   ├── figure_reproduction.py  # Code to reproduce figures from the paper  
│   ├── analysis.py             # Statistical analysis and critical evaluation  
│   ├── utils.py                # Helper functions (e.g., graph generation, energy calculation)  
│   ├── calculate.py            # Core functions for energy and network dynamics  
├── 📂 data                     # Simulation results and precomputed datasets  
├── 📂 notebooks                # Jupyter notebooks for exploratory analysis  
├── 📂 figures                  # Generated plots and comparison with the paper  
├── 📜 requirements.txt         # Dependencies list  
└── 📜 LICENSE                  # License information  
```  

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
python src/monte_carlo_sim.py
```
### **3️⃣ Reproduce Figures from the Paper**  
```bash
python src/figure_reproduction.py
```
### **4️⃣ Explore Results in Jupyter Notebook**  
```bash
jupyter notebook notebooks/analysis.ipynb
```

## **📈 Expected Outcomes & Open Questions**  
🔹 Does **social tension minimization** actually lead to jammed states as described?  
🔹 How does **network structure** (e.g., scale-free vs. complete) affect epidemic outcomes?  
🔹 Can **small rule modifications** (e.g., recovery, delayed reactions) significantly alter conclusions?  
🔹 Are **Monte Carlo update rules (p = 0.5 for neutral moves)** justified?  
🔹 Can we **bridge the gap** between this theoretical model and real-world epidemic data?  

## **🛠 Contributors**  
- **[Your Name]** - Main development & analysis  
- **[Your Collaborators]** - Contributions in simulations, theory, and review  
- **[Supervisor/Professor, if applicable]**  

## **📜 License**  
This project is licensed under the **MIT License** – see the [LICENSE](LICENSE) file for details.  

