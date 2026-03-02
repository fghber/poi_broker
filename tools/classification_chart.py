import numpy as np
from bokeh.plotting import figure, show
from bokeh.models import ColumnDataSource, HoverTool, PolarTransform, LabelSet, Title
from bokeh.colors import RGB
from bokeh.io import output_file

# 1. Define dimensions from your SQL schema
categories = [
    'CV/NOVA', 'E', 'LPV', 'PULS', 
    'PERIODIC(OTHER)', 'QUAS', 'SN', 'YSO'
]
n_vars = len(categories)

# 2. Setup angles (2*pi divided by 8 dimensions)
angles = np.linspace(0, 2*np.pi, n_vars, endpoint=False)

# 3. Sample Data (Simulating rows from your classification table)
# alert_id would be your PK
data = {'alert_id': '2488112286015015008', 'values': [0.85, 0.10, 0.05, 0.02, 0.03, 0.40, 0.10, 0.05]}

max_index = np.argmax(data['values']) 
max_value = data['values'][max_index]
classified_category = categories[max_index]


# Create the figure
p = figure(
    width=700, height=700, 
    title="Classification Probability Radar",
    x_range=(-1.5, 1.5), y_range=(-1.5, 1.5),
    tools="reset,save,wheel_zoom"
)

# add a subtitle
p.add_layout(Title(text="Classified as " + classified_category + "("+str(max_value)+")", text_font_size="10pt", text_font_style="italic"), 'above')

# Use PolarTransform for easy radius/angle mapping
polar_obj = PolarTransform()

# 4. Helper to close the loop (8 points -> 9 points)
angles_closed = np.append(angles, angles[0])
categories_closed = categories + [categories[0]]

# 5. Plot each Alert 
# Append the first value to the end to match the 9-point 'angles_closed'
vals_closed = np.append(data['values'], data['values'][0])
renderers = []
source = ColumnDataSource(data=dict(
    radius=vals_closed,
    angle=angles_closed,
    cat=categories_closed,
    # Ensure this list is also length 9
    alert=[data['alert_id']] * len(vals_closed)
))

color = RGB(181, 137, 0)

# Draw the area
p.patch(x=polar_obj.x, y=polar_obj.y, source=source,
        fill_color=color, fill_alpha=0.15, line_color=color, 
        line_width=2, legend_label=data['alert_id'])

# Draw the points for HoverTool
pt = p.scatter(x=polar_obj.x, y=polar_obj.y, source=source,
                size=10, color=color, alpha=0.6, 
                hover_fill_color="white", legend_label=data['alert_id'])
renderers.append(pt)

# 6. Radial Grid & Axis Labels
radii = [0.2, 0.4, 0.6, 0.8, 1.0]
for r in radii:
    # Circular grid line
    p.circle(0, 0, radius=r, fill_color=None, line_color="grey", line_alpha=0.2)
    # Radial tick labels
    p.text(x=0, y=r, text=[f"{r}"], text_font_size="9pt", 
           text_align="center", text_baseline="middle", text_alpha=0.5)

# 7. Dimension Spokes & Labels
for angle, cat in zip(angles, categories):
    # Draw spoke line
    p.line(x=polar_obj.x, y=polar_obj.y, line_color="grey", line_alpha=0.2,
           source=ColumnDataSource(dict(radius=[0, 1], angle=[angle, angle])))
    
    # Place Category labels slightly outside the 1.0 radius
    label_source = ColumnDataSource(dict(radius=[1.2], angle=[angle], names=[cat]))
    p.text(x=polar_obj.x, y=polar_obj.y, text='names', source=label_source,
           text_align="center", text_baseline="middle", 
           text_font_style="bold", text_font_size="10pt")

# 8. Interactive Features
hover = HoverTool(
    renderers=renderers,
    tooltips=[
        ('Alert ID', '@alert'),
        ('Class', '@cat'),
        ('Probability', '@radius{0.00}')
    ]
)
p.add_tools(hover)

# Styling cleanup
p.grid.grid_line_color = None
p.xaxis.visible = False
p.yaxis.visible = False
#p.legend.click_policy = "hide"
p.legend.location = "top_right"

show(p)