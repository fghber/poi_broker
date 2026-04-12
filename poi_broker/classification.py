from flask import Blueprint, render_template, abort, request, make_response
import numpy as np
from sqlalchemy import text
from bokeh.plotting import figure, show
from bokeh.models import ColumnDataSource, HoverTool, PolarTransform, LabelSet, Title
from bokeh.colors import RGB
from bokeh.embed import components

classification_blueprint = Blueprint('classification', __name__)
# IDEA: Perhaps add a URL prefix classification/

@classification_blueprint.route('/query_classification')
def classification_plot():

    alertId = request.args.get('alertId')
    if not alertId:
        return ('<div class="no-data alert alert-warning">Missing alertId</div>', '') 

    # load values from the SQLite 'classification' table using SQL (no model required)
    # local import to avoid circular import at module import time
    from . import db
    # SQLAlchemy's text() should handle parameterization safely
    sql = text("""
        SELECT p_cvnova, p_e, p_lpv, p_puls,
               p_periodic_other, p_quas, p_sn, p_yso
        FROM classification
        WHERE alert_id = :id
    """)
    row = db.session.execute(sql, {'id': alertId}).fetchone()
    if row is None:
        return (f'<div class="no-data alert alert-warning">No classification data found for alert_id={alertId}</div>', '') 

    # ensure numeric values and replace NULL with 0.0
    values = [float(v) if v is not None else 0.0 for v in row]

    # map to the expected structure used later in the function
    data = {'alert_id': alertId, 'values': values}

    max_index = np.argmax(data['values'])
    max_value = data['values'][max_index]
    classified_category = [
        'CV/NOVA', 'E', 'LPV', 'PULS',
        'PERIODIC(OTHER)', 'QUAS', 'SN', 'YSO'
    ][max_index]


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
    angles = np.linspace(0, 2*np.pi, len(data['values']), endpoint=False)
    angles_closed = np.append(angles, angles[0])
    categories = [
        'CV/NOVA', 'E', 'LPV', 'PULS',
        'PERIODIC(OTHER)', 'QUAS', 'SN', 'YSO'
    ]
    categories_closed = categories + [categories[0]]

    # 5. Plot each Alert 
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

    #show(p)

    # Get Classification Chart Components
    script, div = components(p)

    if not p.renderers:
        return ('<div class="no-data alert alert-warning">There is no classification data to plot!</div>', '')
 
    # Return the components to the HTML template
    return f'{ div }{ script }'
