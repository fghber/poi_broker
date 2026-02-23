from sqlalchemy import inspect

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
    #number with optional deciaml point
    regex = r"[<>]?[+-]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))"
    matches = re.findall(regex, text)
    if len(matches) < 1:
        return None
    elif len(matches) == 1:
        return [matches[0]]
    else:
        return list(map(lambda m: m.replace('>', '').replace('<', ''), matches[0:2]))
         #TODO: we either want to remove </> or add them if missing to be consistent (TBD)

def extract_float_filter(input_field, db_field, query):
    float_func = lambda x: float(x)
    return extract_filter(input_field, db_field, query, float_func)

def extract_int_filter(input_field, db_field, query):
    int_func = lambda x: int(x)
    return extract_filter(input_field, db_field, query, int_func)

def extract_filter(input_field, db_field, query, convert_callback):
    #pdb.set_trace()
    if len(input_field) == 1:
        if '>' in input_field[0]:
            query = query.filter(db_field >= convert_callback(input_field[0].replace('>', '')))
        elif '<' in input_field[0]:
            query = query.filter(db_field <= convert_callback(input_field[0].replace('<', '')))
        else:
            query = query.filter(db_field == convert_callback(input_field[0]))
    else: #2 inputs
        input_field.sort()  #REM: Ensure >min <max order
        print(input_field[0])
        print(input_field[1])
        query = query.filter(db_field >= convert_callback(input_field[0]))
        query = query.filter(db_field <= convert_callback(input_field[1]))
    return query

# Helper function for safe serialization
def safe_serialize(obj, locusId):
    """Safely serialize a dictionary to JSON."""
    try:
        print(dict)
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
