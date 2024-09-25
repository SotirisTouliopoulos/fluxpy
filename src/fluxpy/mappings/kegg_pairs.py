
import re

def extract_dynamic_value(html_string):
    # Extract the dynamic value (e.g., RC00160) using regex
    match = re.search(r'&gt;([^<]+)</a>', html_string)
    if match:
        return match.group(1)
    return None

def transform_values(values):
    if isinstance(values, list) and len(values) == 1 and isinstance(values[0], str):
        dynamic_value = extract_dynamic_value(values[0])
        if dynamic_value:
            return dynamic_value
    return values

def transform_dict(d):
    if isinstance(d, dict):
        transformed_dict = {}
        for key, value in d.items():
            if key == 'values':
                transformed_dict[key] = transform_values(value)
            else:
                transformed_dict[key] = transform_dict(value)
        return transformed_dict
    elif isinstance(d, list):
        return [transform_dict(item) for item in d]
    else:
        return