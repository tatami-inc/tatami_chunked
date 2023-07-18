<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.9.5">
  <compound kind="file">
    <name>CustomChunkedMatrix.hpp</name>
    <path>/github/workspace/include/tatami_chunked/</path>
    <filename>CustomChunkedMatrix_8hpp.html</filename>
    <class kind="struct">tatami_chunked::CustomChunkedOptions</class>
    <class kind="class">tatami_chunked::CustomChunkedDenseMatrix</class>
    <class kind="class">tatami_chunked::CustomChunkedSparseMatrix</class>
    <namespace>tatami_chunked</namespace>
  </compound>
  <compound kind="file">
    <name>LruSlabCache.hpp</name>
    <path>/github/workspace/include/tatami_chunked/</path>
    <filename>LruSlabCache_8hpp.html</filename>
    <class kind="class">tatami_chunked::LruSlabCache</class>
    <namespace>tatami_chunked</namespace>
  </compound>
  <compound kind="file">
    <name>OracleSlabCache.hpp</name>
    <path>/github/workspace/include/tatami_chunked/</path>
    <filename>OracleSlabCache_8hpp.html</filename>
    <class kind="class">tatami_chunked::OracleSlabCache</class>
    <namespace>tatami_chunked</namespace>
  </compound>
  <compound kind="file">
    <name>simple_chunk_wrappers.hpp</name>
    <path>/github/workspace/include/tatami_chunked/</path>
    <filename>simple__chunk__wrappers_8hpp.html</filename>
    <class kind="struct">tatami_chunked::SimpleDenseChunkWrapper</class>
    <class kind="struct">tatami_chunked::SimpleSparseChunkWrapper</class>
    <class kind="struct">tatami_chunked::SimpleSparseChunkWrapper::Workspace</class>
    <namespace>tatami_chunked</namespace>
  </compound>
  <compound kind="file">
    <name>SubsettedOracleSlabCache.hpp</name>
    <path>/github/workspace/include/tatami_chunked/</path>
    <filename>SubsettedOracleSlabCache_8hpp.html</filename>
    <class kind="class">tatami_chunked::SubsettedOracleSlabCache</class>
    <class kind="struct">tatami_chunked::SubsettedOracleSlabCache::SubsetDetails</class>
    <class kind="struct">tatami_chunked::SubsettedOracleSlabCache::CachedSlab</class>
    <namespace>tatami_chunked</namespace>
    <member kind="enumeration">
      <type></type>
      <name>SubsetSelection</name>
      <anchorfile>namespacetatami__chunked.html</anchorfile>
      <anchor>a04acaabc7d3d46a03f1a2faca6dc8ff6</anchor>
      <arglist></arglist>
      <enumvalue file="namespacetatami__chunked.html" anchor="a04acaabc7d3d46a03f1a2faca6dc8ff6aba7de5bc6888294e5884b024a4c894f1">FULL</enumvalue>
      <enumvalue file="namespacetatami__chunked.html" anchor="a04acaabc7d3d46a03f1a2faca6dc8ff6a4d34f53389ed7f28ca91fc31ea360a66">BLOCK</enumvalue>
      <enumvalue file="namespacetatami__chunked.html" anchor="a04acaabc7d3d46a03f1a2faca6dc8ff6acb4ae3b37047fb4b2c0d16f8bf84f076">INDEX</enumvalue>
    </member>
  </compound>
  <compound kind="file">
    <name>tatami_chunked.hpp</name>
    <path>/github/workspace/include/tatami_chunked/</path>
    <filename>tatami__chunked_8hpp.html</filename>
    <namespace>tatami_chunked</namespace>
  </compound>
  <compound kind="file">
    <name>typical_slab_cache.hpp</name>
    <path>/github/workspace/include/tatami_chunked/</path>
    <filename>typical__slab__cache_8hpp.html</filename>
    <class kind="struct">tatami_chunked::TypicalSlabCacheOptions</class>
    <class kind="struct">tatami_chunked::TypicalSlabCacheWorkspace</class>
    <namespace>tatami_chunked</namespace>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::SubsettedOracleSlabCache::CachedSlab</name>
    <filename>structtatami__chunked_1_1SubsettedOracleSlabCache_1_1CachedSlab.html</filename>
    <member kind="variable">
      <type>Slab_</type>
      <name>contents</name>
      <anchorfile>structtatami__chunked_1_1SubsettedOracleSlabCache_1_1CachedSlab.html</anchorfile>
      <anchor>a4c00ff419ab443f384e9ee21139bcc10</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>SubsetDetails</type>
      <name>subset</name>
      <anchorfile>structtatami__chunked_1_1SubsettedOracleSlabCache_1_1CachedSlab.html</anchorfile>
      <anchor>a3ccaae84918e9054d927a630cbf0e006</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami_chunked::CustomChunkedDenseMatrix</name>
    <filename>classtatami__chunked_1_1CustomChunkedDenseMatrix.html</filename>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>typename Chunk_</templarg>
    <templarg>bool use_subsetted_oracle_</templarg>
    <base>VirtualDenseMatrix&lt; Value_, Index_ &gt;</base>
    <member kind="function">
      <type></type>
      <name>CustomChunkedDenseMatrix</name>
      <anchorfile>classtatami__chunked_1_1CustomChunkedDenseMatrix.html</anchorfile>
      <anchor>a105e52c23a1e46681369e0ea48afc591</anchor>
      <arglist>(Index_ mat_nrow, Index_ mat_ncol, Index_ chunk_nrow, Index_ chunk_ncol, std::vector&lt; Chunk_ &gt; chunks, bool row_major, const CustomChunkedOptions &amp;opt)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::CustomChunkedOptions</name>
    <filename>structtatami__chunked_1_1CustomChunkedOptions.html</filename>
    <base>tatami_chunked::TypicalSlabCacheOptions</base>
  </compound>
  <compound kind="class">
    <name>tatami_chunked::CustomChunkedSparseMatrix</name>
    <filename>classtatami__chunked_1_1CustomChunkedSparseMatrix.html</filename>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>typename Chunk_</templarg>
    <templarg>bool use_subsetted_oracle_</templarg>
    <base>Matrix&lt; Value_, Index_ &gt;</base>
    <member kind="function">
      <type></type>
      <name>CustomChunkedSparseMatrix</name>
      <anchorfile>classtatami__chunked_1_1CustomChunkedSparseMatrix.html</anchorfile>
      <anchor>a398ed1cddde3ac16a7ba701a4a4eedc1</anchor>
      <arglist>(Index_ mat_nrow, Index_ mat_ncol, Index_ chunk_nrow, Index_ chunk_ncol, std::vector&lt; Chunk_ &gt; chunks, bool row_major, const CustomChunkedOptions &amp;opt)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami_chunked::LruSlabCache</name>
    <filename>classtatami__chunked_1_1LruSlabCache.html</filename>
    <templarg>typename Id_</templarg>
    <templarg>class Slab_</templarg>
    <member kind="function">
      <type></type>
      <name>LruSlabCache</name>
      <anchorfile>classtatami__chunked_1_1LruSlabCache.html</anchorfile>
      <anchor>ac305164cd47db9f7f008266ddef9f4fb</anchor>
      <arglist>(size_t m)</arglist>
    </member>
    <member kind="function">
      <type>const Slab_ &amp;</type>
      <name>find</name>
      <anchorfile>classtatami__chunked_1_1LruSlabCache.html</anchorfile>
      <anchor>a022825f5749bba03b112fbe72d98ddd5</anchor>
      <arglist>(Id_ id, Cfunction_ create, Pfunction_ populate)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami_chunked::OracleSlabCache</name>
    <filename>classtatami__chunked_1_1OracleSlabCache.html</filename>
    <templarg>typename Id_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>class Slab_</templarg>
    <member kind="function">
      <type></type>
      <name>OracleSlabCache</name>
      <anchorfile>classtatami__chunked_1_1OracleSlabCache.html</anchorfile>
      <anchor>ab2ec1decfa4cdbff71b673d1ef960488</anchor>
      <arglist>(std::unique_ptr&lt; tatami::Oracle&lt; Index_ &gt; &gt; oracle, size_t per_iteration, size_t num_slabs)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; const Slab_ *, Index_ &gt;</type>
      <name>next</name>
      <anchorfile>classtatami__chunked_1_1OracleSlabCache.html</anchorfile>
      <anchor>a8a7fd554334053ebef56bdcb5d6c0321</anchor>
      <arglist>(Ifunction_ identify, Cfunction_ create, Pfunction_ populate)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::SimpleDenseChunkWrapper</name>
    <filename>structtatami__chunked_1_1SimpleDenseChunkWrapper.html</filename>
    <templarg>class SimpleChunk_</templarg>
    <member kind="typedef">
      <type>SimpleChunk_::value_type</type>
      <name>value_type</name>
      <anchorfile>structtatami__chunked_1_1SimpleDenseChunkWrapper.html</anchorfile>
      <anchor>ae5bf9f551bb7b3bd9b8663522b6e4f4c</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>std::vector&lt; value_type &gt;</type>
      <name>Workspace</name>
      <anchorfile>structtatami__chunked_1_1SimpleDenseChunkWrapper.html</anchorfile>
      <anchor>a4c87a7880792f8f92402aec85b57f30c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SimpleDenseChunkWrapper</name>
      <anchorfile>structtatami__chunked_1_1SimpleDenseChunkWrapper.html</anchorfile>
      <anchor>a6849ffc525c47c55d3d76f5ebb894160</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SimpleDenseChunkWrapper</name>
      <anchorfile>structtatami__chunked_1_1SimpleDenseChunkWrapper.html</anchorfile>
      <anchor>ad5b601444cc5a3a8bf1301faa42558a7</anchor>
      <arglist>(SimpleChunk_ c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1SimpleDenseChunkWrapper.html</anchorfile>
      <anchor>ad9ed0914c50251a78a0566f2f0a166ca</anchor>
      <arglist>(Index_ primary, Index_ secondary_start, Index_ secondary_length, Workspace &amp;work, Output_ *output) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1SimpleDenseChunkWrapper.html</anchorfile>
      <anchor>ace0d720292d4b4bad6acd2b219addbf8</anchor>
      <arglist>(Index_ primary, const std::vector&lt; Index_ &gt; &amp;secondary_indices, Workspace &amp;work, Output_ *output) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1SimpleDenseChunkWrapper.html</anchorfile>
      <anchor>a429fbb4e5588ee5a67256222f42a3dd7</anchor>
      <arglist>(Index_ primary_start, Index_ primary_length, Index_ secondary_start, Index_ secondary_length, Workspace &amp;work, Output_ *output, size_t stride) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1SimpleDenseChunkWrapper.html</anchorfile>
      <anchor>a7f06978c01dc2de3e8ba10b42b8c91d6</anchor>
      <arglist>(Index_ primary_start, Index_ primary_length, const std::vector&lt; Index_ &gt; &amp;secondary_indices, Workspace &amp;work, Output_ *output, size_t stride) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1SimpleDenseChunkWrapper.html</anchorfile>
      <anchor>a00faab69eadbc9c2a5850f3f96b4f00f</anchor>
      <arglist>(const std::vector&lt; Index_ &gt; &amp;primary_indices, Index_ secondary_start, Index_ secondary_length, Workspace &amp;work, Output_ *output, size_t stride) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1SimpleDenseChunkWrapper.html</anchorfile>
      <anchor>a86a4ffdad8a474113d04941dab9e201f</anchor>
      <arglist>(const std::vector&lt; Index_ &gt; &amp;primary_indices, const std::vector&lt; Index_ &gt; &amp;secondary_indices, Workspace &amp;work, Output_ *output, size_t stride) const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::SimpleSparseChunkWrapper</name>
    <filename>structtatami__chunked_1_1SimpleSparseChunkWrapper.html</filename>
    <templarg>class SimpleChunk_</templarg>
    <class kind="struct">tatami_chunked::SimpleSparseChunkWrapper::Workspace</class>
    <member kind="typedef">
      <type>SimpleChunk_::value_type</type>
      <name>value_type</name>
      <anchorfile>structtatami__chunked_1_1SimpleSparseChunkWrapper.html</anchorfile>
      <anchor>a8550eae88e491cdc9891de5d5f65f0c9</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>SimpleChunk_::index_type</type>
      <name>index_type</name>
      <anchorfile>structtatami__chunked_1_1SimpleSparseChunkWrapper.html</anchorfile>
      <anchor>aee7f0271689b94d4c76891aa87e4f7a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SimpleSparseChunkWrapper</name>
      <anchorfile>structtatami__chunked_1_1SimpleSparseChunkWrapper.html</anchorfile>
      <anchor>aa1bb067a0363059be4fae4609bf0da40</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SimpleSparseChunkWrapper</name>
      <anchorfile>structtatami__chunked_1_1SimpleSparseChunkWrapper.html</anchorfile>
      <anchor>ad426667f705a2518ffc9199fcf2915d2</anchor>
      <arglist>(SimpleChunk_ c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1SimpleSparseChunkWrapper.html</anchorfile>
      <anchor>a8a062c1bf9cea336bf0dda0e91f6067e</anchor>
      <arglist>(Index_ primary, Index_ secondary_start, Index_ secondary_length, Workspace &amp;work, std::vector&lt; OutputValue_ &gt; &amp;output_values, std::vector&lt; OutputIndex_ &gt; &amp;output_indices, OutputIndex_ shift) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1SimpleSparseChunkWrapper.html</anchorfile>
      <anchor>aa27a8fa69e1c6400d38585f2b0e5a361</anchor>
      <arglist>(Index_ primary, const std::vector&lt; Index_ &gt; &amp;secondary_indices, Workspace &amp;work, std::vector&lt; OutputValue_ &gt; &amp;output_values, std::vector&lt; OutputIndex_ &gt; &amp;output_indices, OutputIndex_ shift) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1SimpleSparseChunkWrapper.html</anchorfile>
      <anchor>a928e327793a4ed9fccc3e18eebdb23c7</anchor>
      <arglist>(Index_ primary_start, Index_ primary_length, Index_ secondary_start, Index_ secondary_length, Workspace &amp;work, std::vector&lt; std::vector&lt; OutputValue_ &gt; &gt; &amp;output_values, std::vector&lt; std::vector&lt; OutputIndex_ &gt; &gt; &amp;output_indices, OutputIndex_ shift) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1SimpleSparseChunkWrapper.html</anchorfile>
      <anchor>a458cd60cd17c81407366041ff2fa9569</anchor>
      <arglist>(Index_ primary_start, Index_ primary_length, const std::vector&lt; Index_ &gt; &amp;secondary_indices, Workspace &amp;work, std::vector&lt; std::vector&lt; OutputValue_ &gt; &gt; &amp;output_values, std::vector&lt; std::vector&lt; OutputIndex_ &gt; &gt; &amp;output_indices, OutputIndex_ shift) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1SimpleSparseChunkWrapper.html</anchorfile>
      <anchor>a21aeb9076df2cb030520052b14b5bae7</anchor>
      <arglist>(const std::vector&lt; Index_ &gt; &amp;primary_indices, Index_ secondary_start, Index_ secondary_length, Workspace &amp;work, std::vector&lt; std::vector&lt; OutputValue_ &gt; &gt; &amp;output_values, std::vector&lt; std::vector&lt; OutputIndex_ &gt; &gt; &amp;output_indices, OutputIndex_ shift) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>extract</name>
      <anchorfile>structtatami__chunked_1_1SimpleSparseChunkWrapper.html</anchorfile>
      <anchor>a5c746e76474d14e55c7bcdabd06dfd50</anchor>
      <arglist>(const std::vector&lt; Index_ &gt; &amp;primary_indices, const std::vector&lt; Index_ &gt; &amp;secondary_indices, Workspace &amp;work, std::vector&lt; std::vector&lt; OutputValue_ &gt; &gt; &amp;output_values, std::vector&lt; std::vector&lt; OutputIndex_ &gt; &gt; &amp;output_indices, OutputIndex_ shift) const</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::SubsettedOracleSlabCache::SubsetDetails</name>
    <filename>structtatami__chunked_1_1SubsettedOracleSlabCache_1_1SubsetDetails.html</filename>
    <member kind="variable">
      <type>SubsetSelection</type>
      <name>selection</name>
      <anchorfile>structtatami__chunked_1_1SubsettedOracleSlabCache_1_1SubsetDetails.html</anchorfile>
      <anchor>a65b521f2d14bc202f412da9d7f9ca92f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Index_</type>
      <name>block_start</name>
      <anchorfile>structtatami__chunked_1_1SubsettedOracleSlabCache_1_1SubsetDetails.html</anchorfile>
      <anchor>a8a338261ab5d0763da698ebc295704cd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Index_</type>
      <name>block_length</name>
      <anchorfile>structtatami__chunked_1_1SubsettedOracleSlabCache_1_1SubsetDetails.html</anchorfile>
      <anchor>a98d7097252e23c52e81f07330dde84e6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Index_ &gt;</type>
      <name>indices</name>
      <anchorfile>structtatami__chunked_1_1SubsettedOracleSlabCache_1_1SubsetDetails.html</anchorfile>
      <anchor>a9f8cd50f3419f1620b7765ec56ab9c27</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::unordered_map&lt; Index_, Index_ &gt;</type>
      <name>mapping</name>
      <anchorfile>structtatami__chunked_1_1SubsettedOracleSlabCache_1_1SubsetDetails.html</anchorfile>
      <anchor>a92222ed5250f12b0fcb58b7d7251658b</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>tatami_chunked::SubsettedOracleSlabCache</name>
    <filename>classtatami__chunked_1_1SubsettedOracleSlabCache.html</filename>
    <templarg>typename Id_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>class Slab_</templarg>
    <class kind="struct">tatami_chunked::SubsettedOracleSlabCache::CachedSlab</class>
    <class kind="struct">tatami_chunked::SubsettedOracleSlabCache::SubsetDetails</class>
    <member kind="function">
      <type></type>
      <name>SubsettedOracleSlabCache</name>
      <anchorfile>classtatami__chunked_1_1SubsettedOracleSlabCache.html</anchorfile>
      <anchor>a4d007072fdee018b8cd075784759dc19</anchor>
      <arglist>(std::unique_ptr&lt; tatami::Oracle&lt; Index_ &gt; &gt; oracle, size_t per_iteration, size_t num_slabs)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; const CachedSlab *, Index_ &gt;</type>
      <name>next</name>
      <anchorfile>classtatami__chunked_1_1SubsettedOracleSlabCache.html</anchorfile>
      <anchor>ac64bbad23fe18047559614afac774728</anchor>
      <arglist>(Ifunction_ identify, Cfunction_ create, Pfunction_ populate)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::TypicalSlabCacheOptions</name>
    <filename>structtatami__chunked_1_1TypicalSlabCacheOptions.html</filename>
    <member kind="variable">
      <type>size_t</type>
      <name>maximum_cache_size</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheOptions.html</anchorfile>
      <anchor>acb6a39320c2eccd7cf2231e8b29b13fc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>require_minimum_cache</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheOptions.html</anchorfile>
      <anchor>ac48f95c0943d497e3d743b9573d47510</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::TypicalSlabCacheWorkspace</name>
    <filename>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</filename>
    <templarg>typename Index_</templarg>
    <templarg>class Slab_</templarg>
    <templarg>bool subset_</templarg>
    <member kind="typedef">
      <type>std::conditional&lt; subset_, SubsettedOracleSlabCache&lt; Index_, Index_, Slab_ &gt;, OracleSlabCache&lt; Index_, Index_, Slab_ &gt; &gt;::type</type>
      <name>OracleCache</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>a6a4fa7d2a4fc66947a3646ea9a50948c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TypicalSlabCacheWorkspace</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>ac27a71aa118a875f1f915e7dfc0dfbd7</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TypicalSlabCacheWorkspace</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>a347e990db1a0471ea42a1f27702d4117</anchor>
      <arglist>(Index_ primary_length, Index_ secondary_length, size_t cache_size_in_elements, bool require_minimum_cache)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_oracle</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>a44ac93763f6c372657867e410f0dc80b</anchor>
      <arglist>(std::unique_ptr&lt; tatami::Oracle&lt; Index_ &gt; &gt; o)</arglist>
    </member>
    <member kind="variable">
      <type>Index_</type>
      <name>primary_length</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>aa1f492cff84927346fe0c2da318f1e0c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>size_t</type>
      <name>slab_size_in_elements</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>a95f2fbaa94e3c78264b4ad6b742e3dd0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>size_t</type>
      <name>num_slabs_in_cache</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>a9b01d451085734c970f4c6f985443cbf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::unique_ptr&lt; LruSlabCache&lt; Index_, Slab_ &gt; &gt;</type>
      <name>lru_cache</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>af21e3f5737b73c1155ef12a7a9896269</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::unique_ptr&lt; OracleCache &gt;</type>
      <name>oracle_cache</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>a9dc2bab8876b8d7be33337de133ca8a7</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami_chunked::SimpleSparseChunkWrapper::Workspace</name>
    <filename>structtatami__chunked_1_1SimpleSparseChunkWrapper_1_1Workspace.html</filename>
  </compound>
  <compound kind="namespace">
    <name>tatami_chunked</name>
    <filename>namespacetatami__chunked.html</filename>
    <class kind="class">tatami_chunked::CustomChunkedDenseMatrix</class>
    <class kind="struct">tatami_chunked::CustomChunkedOptions</class>
    <class kind="class">tatami_chunked::CustomChunkedSparseMatrix</class>
    <class kind="class">tatami_chunked::LruSlabCache</class>
    <class kind="class">tatami_chunked::OracleSlabCache</class>
    <class kind="struct">tatami_chunked::SimpleDenseChunkWrapper</class>
    <class kind="struct">tatami_chunked::SimpleSparseChunkWrapper</class>
    <class kind="class">tatami_chunked::SubsettedOracleSlabCache</class>
    <class kind="struct">tatami_chunked::TypicalSlabCacheOptions</class>
    <class kind="struct">tatami_chunked::TypicalSlabCacheWorkspace</class>
    <member kind="enumeration">
      <type></type>
      <name>SubsetSelection</name>
      <anchorfile>namespacetatami__chunked.html</anchorfile>
      <anchor>a04acaabc7d3d46a03f1a2faca6dc8ff6</anchor>
      <arglist></arglist>
      <enumvalue file="namespacetatami__chunked.html" anchor="a04acaabc7d3d46a03f1a2faca6dc8ff6aba7de5bc6888294e5884b024a4c894f1">FULL</enumvalue>
      <enumvalue file="namespacetatami__chunked.html" anchor="a04acaabc7d3d46a03f1a2faca6dc8ff6a4d34f53389ed7f28ca91fc31ea360a66">BLOCK</enumvalue>
      <enumvalue file="namespacetatami__chunked.html" anchor="a04acaabc7d3d46a03f1a2faca6dc8ff6acb4ae3b37047fb4b2c0d16f8bf84f076">INDEX</enumvalue>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>tatami bindings for chunked matrices</title>
    <filename>index.html</filename>
    <docanchor file="index.html">md__github_workspace_README</docanchor>
  </compound>
</tagfile>
