# PKE Documentation

[ciphertext.h](ciphertext.h)

- for the representation of ciphertext in OpenFHE

- provides `CiphertextImpl` which is used to contain encrypted text

[ciphertext-ser.h](ciphertext-ser.h)

- exposes serialization methods for ciphertexts to [USCiLab - cereal](https://github.com/USCiLab/cereal)

- must be included any time we need ciphertext serialization

[constants.h](constants.h)

- Contains the various constants used throughout the `PKE` module including:
  - PKE Scheme Feature: Lists all features supported by PKE schemees.
    - PKE (PKE)
    - Keyswitch (KEYSWITCH)
    - proxy-reencryption (PRE)
    - Leveled Somewhat Homomorphic Encryption (LEVELED SHE)
    - Advanced Somewhat Homomorphic Encryption (ADVANCEDSHE)
    - Threshold FHE (MULTIPARTY)
    - Fully Homomorphic Encryption (FHE)
  
  - MODE 
    - Ring Learning with Error (RLWE)
    - Optimized (OPTIMIZED)
    - Sparse (sparse)
    
  - Rescaling Technique 
    - Fixed Manual (FIXEDMANUAL)
    - Fixed Auto (FIXEDAUTO)
    - Flexible Auto (FLEXIBLEAUTO)
    - No Rescaling (NORESCALE)

  - Key Switch Technique 
    - BV (See [Fully Homomorphic Encryption from Ring-LWE and Security for Key Dependent Messages (BV Scheme)](https://www.wisdom.weizmann.ac.il/~zvikab/localpapers/IdealHom.pdf))
    - Hybrid (See [Homomorphic Evaluation of the AES Circuit](https://eprint.iacr.org/2012/099.pdf))
    
  - Encryption Technique
    - Standard (STANDARD)
    - P Over Q (POVERQ)
    
  - Multiplication Technique
    - Bajard-Eynard-Hasan-Zucca (BEHZ)
    - Halevi-Polyakov-Shoup (HPS)
    - Halevi-Polyakov-Shoup P over Q (HPSPOVERQ)
    - Halevi-Polyakov-Shoup P over Q Leveled Multiplication (HPSPOVERQLEVELED)

[cryptocontextfactory.h](cryptocontextfactory.h)

- Generates new `CryptoContexts` from user parameters

[cryptocontextgen.h](cryptocontextgen.h)

- `Cryptocontext` generator for the various PKE schemes (notably `BFV`, `BGV`, and `CKKS` in their `RNS` forms)

[cryptocontext.h](cryptocontext.h)

- defines `CryptoContextImpl`, which is used to access the OpenFHE library

- all OpenFHE objects are created "within" a `CryptoContext` which acts like an object "manager". 
Objects can only be used in the context they were created in

[cryptocontexthelper.h](cryptocontexthelper.h)

- provides helper functions to print out the parameter sets 

[cryptocontextparametersets.h](cryptocontextparametersets.h)

- defines the parameter sets for `CryptoContexts`

[cryptocontext-ser.h](cryptocontext-ser.h)

- exposes serialization methods for `CryptoContext` to [USCiLab - cereal](https://github.com/USCiLab/cereal)

- must be included any time we need `CryptoContext` serialization

[cryptoobject.h](cryptoobject.h)

- comprises a `context`, and `keytag`:
  - `context`: `CryptoContext` this object belongs to
  - `keytag`: tag that is used to find the evaluation key needed for various operations

[gen-cryptocontext.h](gen-cryptocontext.h)

- Constructs `CryptoContext` based on the provided set of parameters

[globals.h](globals.h)

- Global value definitions

[metadata.h](metadata.h)

- metadata container and helper function definition

[metadata-ser.h](metadata-ser.h)

- exposes serialization methods for `Metadata` to [USCiLab - cereal](https://github.com/USCiLab/cereal)

- must be included any time we need `Metadata` serialization

[openfhe.h](openfhe.h)

- top level for ease of import