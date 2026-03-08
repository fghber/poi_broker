from sqlalchemy import inspect
import re
import json
from astropy.time import Time

# Helper for converting SQLAlchemy objects to dictionaries
def object_as_dict(obj):
    return {c.key: getattr(obj, c.key)
            for c in inspect(obj).mapper.column_attrs}

def result_to_dict(query_results):
    def to_dict(obj):
        return {c.name: getattr(obj, c.name) for c in obj.__table__.columns}    
    return [to_dict(result) for result in query_results]


# Helper for filtering
def extract_numbers(text):
    #number with optional decimal point
    regex = r"[<>]?[+-]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))"
    matches = re.findall(regex, text)
    if len(matches) < 1:
        return None
    elif len(matches) == 1:
        return [matches[0]]
    else:
        return list(map(lambda m: m.replace('>', '').replace('<', ''), matches[0:2]))
         #TODO: we either want to remove </> or add them if missing to be consistent (TBD)

def extract_dates(text) -> list[str]:
    #date in yyyymmdd format, with optional > or < for filtering
    regex = r"[<>]?\d{8}"
    matches = re.findall(regex, text)
    if len(matches) < 1:
        regex_iso = r"[<>]?\d{4}-\d{2}-\d{2}(?:[ T]\d{2}:\d{2}:\d{2})?"
        matches = re.findall(regex_iso, text)

    if len(matches) < 1:
        return []
    elif len(matches) == 1:
        return [matches[0]]
    else:
        return list(map(lambda m: m.replace('>', '').replace('<', ''), matches[0:2]))

def extract_float_filter(input_field, db_field, query):
    float_func = lambda x: float(x)
    return extract_filter(input_field, db_field, query, float_func)

def extract_int_filter(input_field, db_field, query):
    int_func = lambda x: int(x)
    return extract_filter(input_field, db_field, query, int_func)

def extract_mjd_filter(input_field, db_field, query):
    mjd_func = lambda x: Time(x, format='iso', scale='utc').mjd
    return extract_filter(input_field, db_field, query, mjd_func, offset_mjd=1.0/(3600*24)) # 1 second

def extract_filter(input_field, db_field, query, convert_callback, offset_mjd = 0.0):
    #pdb.set_trace()
    if len(input_field) == 1:
        if '>' in input_field[0]:
            query = query.filter(db_field >= convert_callback(input_field[0].replace('>', '')) - offset_mjd)
        elif '<' in input_field[0]:
            query = query.filter(db_field <= convert_callback(input_field[0].replace('<', '')) + offset_mjd)
        else:
            lower_bound = convert_callback(input_field[0]) - offset_mjd
            upper_bound = convert_callback(input_field[0]) + offset_mjd
            query = query.filter(db_field >= lower_bound)
            query = query.filter(db_field <= upper_bound)
            # MJD time can suffer from rounding errors: DB: 61056.12116899993 vs calc: 61056.12116898148
    else: #2 inputs
        filter_fields = [convert_callback(x) for x in input_field]
        filter_fields.sort()  #NOTE: Ensure >min <max order is correct regardless of input order
        lower_bound = filter_fields[0] - offset_mjd
        upper_bound = filter_fields[1] + offset_mjd
        query = query.filter(db_field >= lower_bound)
        query = query.filter(db_field <= upper_bound)
    return query

# Helper function for safe serialization
def safe_serialize(obj):
    """Safely serialize a dictionary to JSON."""
    try:
        return json.dumps(obj)
    except TypeError as e:
        print(f"Serialization error: {e}")
        return json.dumps(serialize_fallback(obj))
   
def serialize_fallback(obj):
    """Fallback handler to make the object serializable by converting binary data to string."""
    if isinstance(obj, bytes):
        return obj.decode('utf-8')  # Convert binary to string
    elif isinstance(obj, dict):
        return {k: serialize_fallback(v) for k, v in obj.items()}  # Keep dict as-is and process its values
    elif isinstance(obj, list):
        return [serialize_fallback(v) for v in obj]  # Keep list as-is and process its elements
    else:
        return obj  # Return other types as-is
