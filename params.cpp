#include "params.h"

void initialize_params(struct params &_params) {
    if (_params.isInit) return;

    _params.out_vcf.open(_params.output_file);
    _params.isInit = true;
}

void deinitialize_params(struct params &_params) {
    if (!_params.isInit) return;

    _params.out_vcf.close();
    _params.isInit = false;
}
