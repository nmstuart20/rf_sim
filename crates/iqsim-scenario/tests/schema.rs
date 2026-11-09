//! The checked-in JSON Schema matches the Rust types.

#![forbid(unsafe_code)]

use std::path::Path;

#[test]
fn checked_in_schema_is_current() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../schema/scenario.schema.json");
    let checked_in = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("{}: {e}; run `cargo xtask schema`", path.display()));
    assert!(
        checked_in == iqsim_scenario::scenario_schema_json(),
        "schema/scenario.schema.json is stale; run `cargo xtask schema` and commit the result"
    );
}

#[test]
fn schema_is_strict_and_documents_units() {
    let schema: serde_json::Value =
        serde_json::from_str(&iqsim_scenario::scenario_schema_json()).unwrap();
    assert_eq!(schema["additionalProperties"], false);
    let defs = &schema["$defs"];
    for strict in [
        "Platform", "Emitter", "Receiver", "Channel", "Sink", "Resource", "Metadata",
    ] {
        assert_eq!(defs[strict]["additionalProperties"], false, "{strict}");
    }
    assert!(
        defs["Frequency"]["pattern"]
            .as_str()
            .unwrap()
            .contains("MHz")
    );
    assert!(defs["Position"]["oneOf"].is_array());
    let text = iqsim_scenario::scenario_schema_json();
    assert!(
        !text.contains("null"),
        "TOML has no null; optional keys are omitted"
    );
    assert_eq!(
        defs["Emitter"]["properties"]["mount"]["$ref"],
        "#/$defs/Ident"
    );
    let required = schema["required"].as_array().unwrap();
    for key in ["schema_version", "metadata", "epoch", "seed", "run_mode"] {
        assert!(required.iter().any(|v| v == key), "{key}");
    }
}
