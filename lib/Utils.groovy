class Utils {
    // Reads params.X once into a real boolean. Nextflow keeps returning the
    // original CLI string on every later read of params.X (the CLI value
    // always wins over an in-script reassignment), so this must be called
    // once per value used, not relied on via a global params mutation.
    static boolean asBool(value) {
        return value instanceof String ? value.toBoolean() : (value as boolean)
    }
}
