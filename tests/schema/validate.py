#!/usr/bin/env python3
"""validate.py: check a JSON file against a JSON Schema, standard library only.

  python3 tests/schema/validate.py SCHEMA.json DOC.json

Supports the keywords tests/schema/summary.schema.json uses: type (a name or a
list), enum, required, properties, additionalProperties (a schema), items,
minItems, maxItems, minimum, allOf and local $ref ("#/definitions/..."). It
exits 1 and prints every error when the document does not match, and exits 2
when the schema uses a keyword it does not know (so a check never passes
because a rule was silently ignored).
"""
import json
import sys

KNOWN = {"$schema", "title", "description", "definitions", "type", "enum", "required", "properties",
         "additionalProperties", "items", "minItems", "maxItems", "minimum", "allOf", "$ref"}
TYPES = {
    "object": lambda v: isinstance(v, dict),
    "array": lambda v: isinstance(v, list),
    "string": lambda v: isinstance(v, str),
    "integer": lambda v: isinstance(v, int) and not isinstance(v, bool),
    "number": lambda v: isinstance(v, (int, float)) and not isinstance(v, bool),
    "boolean": lambda v: isinstance(v, bool),
    "null": lambda v: v is None,
}


class SchemaError(Exception):
    pass


def resolve(root, ref):
    if not ref.startswith("#/"):
        raise SchemaError(f"only local $ref is supported: {ref}")
    node = root
    for part in ref[2:].split("/"):
        node = node[part]
    return node


def check(root, schema, value, path, errors):
    unknown = set(schema) - KNOWN
    if unknown:
        raise SchemaError(f"{path}: unsupported keyword(s) {sorted(unknown)}")
    if "$ref" in schema:
        check(root, resolve(root, schema["$ref"]), value, path, errors)
    for sub in schema.get("allOf", []):
        check(root, sub, value, path, errors)
    if "type" in schema:
        types = schema["type"] if isinstance(schema["type"], list) else [schema["type"]]
        if not any(TYPES[t](value) for t in types):
            errors.append(f"{path}: expected {'/'.join(types)}, got {type(value).__name__} {value!r:.60}")
            return
    if "enum" in schema and value not in schema["enum"]:
        errors.append(f"{path}: {value!r:.60} is not one of {schema['enum']}")
    if "minimum" in schema and TYPES["number"](value) and value < schema["minimum"]:
        errors.append(f"{path}: {value} < minimum {schema['minimum']}")
    if isinstance(value, dict):
        for key in schema.get("required", []):
            if key not in value:
                errors.append(f"{path}: missing required key '{key}'")
        props = schema.get("properties", {})
        for key, sub in props.items():
            if key in value:
                check(root, sub, value[key], f"{path}.{key}", errors)
        extra = schema.get("additionalProperties")
        if isinstance(extra, dict):
            for key in value:
                if key not in props:
                    check(root, extra, value[key], f"{path}.{key}", errors)
    if isinstance(value, list):
        if "minItems" in schema and len(value) < schema["minItems"]:
            errors.append(f"{path}: {len(value)} items < minItems {schema['minItems']}")
        if "maxItems" in schema and len(value) > schema["maxItems"]:
            errors.append(f"{path}: {len(value)} items > maxItems {schema['maxItems']}")
        if "items" in schema:
            for i, item in enumerate(value):
                check(root, schema["items"], item, f"{path}[{i}]", errors)


def validate(schema, doc):
    errors = []
    check(schema, schema, doc, "$", errors)
    return errors


def main(argv):
    if len(argv) != 3:
        print(__doc__.strip().split("\n\n")[1], file=sys.stderr)
        return 2
    with open(argv[1]) as f:
        schema = json.load(f)
    with open(argv[2]) as f:
        doc = json.load(f)
    try:
        errors = validate(schema, doc)
    except SchemaError as e:
        print(f"SCHEMA ERROR: {e}", file=sys.stderr)
        return 2
    for e in errors:
        print(f"INVALID {e}")
    if errors:
        print(f"{argv[2]}: {len(errors)} schema error(s)")
        return 1
    print(f"{argv[2]}: valid against {argv[1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
