import {run_lsoda_tests} from "./lsoda_test";
import {run_lsodar_tests} from "./lsodar_test";

run_all_tests();

function run_all_tests() {
    run_lsoda_tests();
    run_lsodar_tests();
}