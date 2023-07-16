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
  <compound kind="class">
    <name>tatami_chunked::CustomChunkedDenseMatrix</name>
    <filename>classtatami__chunked_1_1CustomChunkedDenseMatrix.html</filename>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>typename Chunk_</templarg>
    <base>VirtualDenseMatrix&lt; Value_, Index_ &gt;</base>
    <member kind="function">
      <type></type>
      <name>CustomChunkedDenseMatrix</name>
      <anchorfile>classtatami__chunked_1_1CustomChunkedDenseMatrix.html</anchorfile>
      <anchor>a22c435aa278294be2af18f81878cea51</anchor>
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
    <base>Matrix&lt; Value_, Index_ &gt;</base>
    <member kind="function">
      <type></type>
      <name>CustomChunkedSparseMatrix</name>
      <anchorfile>classtatami__chunked_1_1CustomChunkedSparseMatrix.html</anchorfile>
      <anchor>ade92e3332c338cfb65f0ba366d2eaa14</anchor>
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
    <member kind="function">
      <type></type>
      <name>TypicalSlabCacheWorkspace</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>ab99537bf021828325a9c9d78cf6f9458</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TypicalSlabCacheWorkspace</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>a116dd3a78cd0ab9f43e99425f43bfa69</anchor>
      <arglist>(Index_ primary_length, Index_ secondary_length, size_t cache_size_in_elements, bool require_minimum_cache)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_oracle</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>a64e99c4866ed080bfbf7dd63ac112474</anchor>
      <arglist>(std::unique_ptr&lt; tatami::Oracle&lt; Index_ &gt; &gt; o)</arglist>
    </member>
    <member kind="variable">
      <type>Index_</type>
      <name>primary_length</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>a46f8deb0af915ddaa539da64efa0010d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>size_t</type>
      <name>slab_size_in_elements</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>a799cf52c298fefc718bf8a8137139933</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>size_t</type>
      <name>num_slabs_in_cache</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>a4db6804025caf7a08ffcce71ac9c36d6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::unique_ptr&lt; LruSlabCache&lt; Index_, Slab_ &gt; &gt;</type>
      <name>lru_cache</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>a064f91ff40a9b14db716b17098061515</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::unique_ptr&lt; OracleSlabCache&lt; Index_, Index_, Slab_ &gt; &gt;</type>
      <name>oracle_cache</name>
      <anchorfile>structtatami__chunked_1_1TypicalSlabCacheWorkspace.html</anchorfile>
      <anchor>a9edf66abbcc4d3c07cd3bee543da2598</anchor>
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
    <class kind="struct">tatami_chunked::TypicalSlabCacheOptions</class>
    <class kind="struct">tatami_chunked::TypicalSlabCacheWorkspace</class>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>tatami bindings for chunked matrices</title>
    <filename>index.html</filename>
    <docanchor file="index.html">md__github_workspace_README</docanchor>
  </compound>
</tagfile>
