


# CoBraPython
Computational fluid dynamics simulator to calculate the thermal performance of a full cooling branch

---

### To install dependencies:
```bash
(sudo) pip3 install -r requirements.txt
```

---

### To specify new conditions:
1. Navigate to the root folder (`CoBraPython`)
2. Create a new xml template: `cp CoBraV1a_example.xml <new_setup>.xml`
3. Modify `<new_setup>.xml` as fit (see [xml parameters](#xml-parameters))
4. Navigate to `FluidDynamics/CoolingBranch_v1a.py` and modify the Mass Flow on line 459 as fit.

---

### To run:
```bash
python3 CoolingBranch_v1a.py <path/to/xml>
```

---

### XML parameters:

**Params:**

| param | type | unit | description |
| - | - | - | - |
| <p align="center">Section</p> | <p align="center">integer</p> | <p align="center">dimensionless </p> | <p align="center">Branch section number (0..N)</p>
| <p align="center">Hydraulic Diameter</p> | <p align="center">float</p> | <p align="center">mm </p> | <p align="center">——</p>
| <p align="center">Tube Section Length</p> | <p align="center">float</p> | <p align="center">m</p> | <p align="center">——</p>
| <p align="center">Tube Orientation Angle</p> | <p align="center">float</p> | <p align="center">degrees</p> | <p align="center">——</p>
| <p align="center">Tube Roughness</p> | <p align="center">float</p> | <p align="center">&mu;m </p> | <p align="center">——</p>
| <p align="center">HX Node #</p> | integer  | <p align="center">dimensionless </p> | <p align="center">This field corresponds to the section that each node is in contact with. If a node is not in contact with any other section, this field should be the section number of the node</p>
| <p align="center">HX Flow Direction</p> | <p align="center">integer </p>  | <p align="center">dimensionless </p> | <p align="center">This field corresponds to the flow direction of the sections in contact. Use `1` for co-current, `-1` for counter-current, `0` if not in contact.</p>
| <p align="center">HX Conductance</p> | <p align="center">float</p>  | <p align="center">![equation](https://latex.codecogs.com/gif.latex?%5Cfrac%7BWatts%7D%7Bmeters*Kelvin%7D)</p> | <p align="center">This field specifies the conductance  between two sections in contact. Put `0` if not in contact.</p>
| <p align="center">HX Thickness</p> | <p align="center">float </p>  | <p align="center">mm </p> | <p align="center">This field corresponds to the thickness of the barrier between two sections in contact. Put `0` if not in contact.</p>
| <p align="center">Insulation Conductance</p> | <p align="center">float </p>  | <p align="center">![equation](https://latex.codecogs.com/gif.latex?%5Cfrac%7BWatts%7D%7Bmeters*Kelvin%7D)</p> | <p align="center">——</p>
| <p align="center">Environmental Temperature</p> | <p align="center">float </p>  | <p align="center">Celsius</p> | <p align="center">——</p>
| <p align="center">Tube Wall Thickness</p> | <p align="center">float </p>  | <p align="center">mm</p> | <p align="center">——</p>
| <p align="center">Insulation Thickness</p> | <p align="center">float </p>  | <p align="center">mm</p> | <p align="center">——</p>
| <p align="center">Tube Thermal Conductance</p> | <p align="center">float </p>  | <p align="center">![equation](https://latex.codecogs.com/gif.latex?%5Cfrac%7BWatts%7D%7Bmeters*Kelvin%7D)</p> | <p align="center">——</p>
| <p align="center">Heat Transfer Coefficient (Air)</p> | <p align="center">float </p>  | <p align="center">![equation](https://latex.codecogs.com/gif.latex?%5Cfrac%7BWatts%7D%7Bmeters%5E2*Kelvin%7D)</p> | <p align="center">——</p>
| <p align="center">Heat Flow</p> | <p align="center">float </p>  | <p align="center">Watts</p> | <p align="center">——</p>
| <p align="center">&Delta;L</p> | <p align="center">float </p>  | <p align="center">m</p> | <p align="center">This field corresponds to the length of each chopped section. Default: ![equation](https://latex.codecogs.com/gif.latex?10%5E%7B-3%7D%5Ctext%7B%20m%7D)</p>

