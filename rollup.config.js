import typescript from 'rollup-plugin-typescript2';
import { terser } from "rollup-plugin-terser";

export default {
    input: './src/odepack.ts',
    output: [
        {
            format: 'umd',
            file: `dist/odepack.min.js`,
            name: 'odepack',
        },
    ],
	plugins: [
        terser(),
		typescript(/*{ plugin options }*/)
    ],
    // external: [
        // "tslib"
    // ]
}