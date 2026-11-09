//! Scenario schema, validation, and resolution.
//! Per the plan section "Interfaces" (Scenario files), this crate owns strict
//! unit-aware configuration and named references to platforms, emitters, receivers,
//! libraries, output sinks, and shared resource capacities.
//!
//! Dependencies: iqsim-core, iqsim-units, iqsim-contracts.
//!
//! Scenario files are TOML. [`parse`] reads one into a [`Scenario`] (every table
//! rejects unknown keys, and every quantity needs a unit), [`validate`] runs the
//! semantic checks, and [`lint`] does both. Diagnostics carry byte spans that
//! [`Diagnostic::render`] turns into `file:line:column` messages.
//! [`scenario_schema_json`] is the JSON Schema checked in at
//! `schema/scenario.schema.json` for editors and adapters.

#![forbid(unsafe_code)]

mod diagnostic;
mod epoch;
mod ident;
mod model;
mod validate;

pub use diagnostic::{Diagnostic, line_col};
pub use epoch::Epoch;
pub use ident::{Ident, LibraryRef};
pub use model::*;
pub use validate::validate;

/// Parses scenario TOML. The first syntax, type or unit error stops parsing.
pub fn parse(source: &str) -> Result<Scenario, Diagnostic> {
    toml::from_str(source).map_err(|e| Diagnostic::from_toml(&e))
}

/// Parses and validates; returns every diagnostic (empty when the scenario is clean).
pub fn lint(source: &str) -> Vec<Diagnostic> {
    match parse(source) {
        Ok(scenario) => validate(&scenario),
        Err(diagnostic) => vec![diagnostic],
    }
}

/// The scenario JSON Schema (draft 2020-12).
///
/// TOML has no null, so the `null` alternatives schemars emits for optional fields
/// are removed: an optional key is simply left out.
pub fn scenario_schema() -> schemars::Schema {
    let mut value = schemars::schema_for!(Scenario).to_value();
    strip_null(&mut value);
    schemars::Schema::try_from(value).expect("still a schema object")
}

/// Rewrites `anyOf: [X, {type: null}]` to `X` and drops `"null"` from type arrays.
fn strip_null(value: &mut serde_json::Value) {
    use serde_json::Value;
    match value {
        Value::Object(map) => {
            let is_null = |v: &Value| v.get("type").is_some_and(|t| t == "null");
            if let Some(Value::Array(options)) = map.get("anyOf")
                && options.len() == 2
                && options.iter().any(is_null)
            {
                let keep = options.iter().find(|v| !is_null(v)).cloned();
                map.remove("anyOf");
                if let Some(Value::Object(keep)) = keep {
                    map.extend(keep);
                }
            }
            if let Some(Value::Array(types)) = map.get_mut("type") {
                types.retain(|t| t != "null");
                if types.len() == 1 {
                    let only = types.remove(0);
                    map.insert("type".to_owned(), only);
                }
            }
            map.values_mut().for_each(strip_null);
        }
        Value::Array(items) => items.iter_mut().for_each(strip_null),
        _ => {}
    }
}

/// The scenario JSON Schema as pretty-printed JSON with a trailing newline, exactly
/// as checked in at `schema/scenario.schema.json`.
pub fn scenario_schema_json() -> String {
    let mut text = serde_json::to_string_pretty(&scenario_schema()).expect("a schema serializes");
    text.push('\n');
    text
}
